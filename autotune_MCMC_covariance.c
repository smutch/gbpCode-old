#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void autotune_MCMC_covariance(MCMC_info *MCMC){
  int            n_P;
  int            n_P_squared;
  double        *P_chain;
  double         threshold;
  double         success;
  int            n_autotune;
  int            n_covariance;
  int            n_success;
  int            i_autotune;
  int            i_iteration;
  int            flag_success;
  double         covariance_threshold;
  double        *P_i_bar_accum;
  double        *P_ij_bar_accum;
  double        *V_compute;
  double        *V_compute_last;
  double         difference;
  double         difference_i;
  int            i_P,j_P,k_P;  
  time_t         time_start, time_stop;
  double         time_diff;

  SID_log("Autotuning the covariance matrix (threshold=%6.2lf%%)...",SID_LOG_OPEN|SID_LOG_TIMER,MCMC->covariance_threshold);

  // Grab some things from the MCMC structure
  n_P    =MCMC->n_P;
  P_chain=MCMC->P_chain;

  // Initialize covariance matrix
  n_P_squared   =n_P*n_P;
  P_i_bar_accum =(double *)SID_malloc(sizeof(double)*n_P);
  P_ij_bar_accum=(double *)SID_malloc(sizeof(double)*n_P_squared);
  V_compute     =(double *)SID_malloc(sizeof(double)*n_P_squared);
  V_compute_last=(double *)SID_malloc(sizeof(double)*n_P_squared);
  for(i_P=0,k_P=0;i_P<n_P;i_P++){
    P_i_bar_accum[i_P]=0.;
    for(j_P=0;j_P<n_P;j_P++,k_P++)
      P_ij_bar_accum[k_P]=0.;
  }
  memcpy(V_compute_last,MCMC->V,sizeof(double)*n_P_squared);

  // Initialize everything else
  covariance_threshold=MCMC->covariance_threshold;
  n_autotune          =MCMC->n_autotune_covariance;
  i_autotune          =0;
  difference          =2.*covariance_threshold;
  n_success           =0;
  n_covariance        =0;

  // Start from the initial conditions for each iteration
  MCMC->flag_init_chain=TRUE;

  // Integrate the covariance matrix until it convergence
  while(difference>covariance_threshold){

    // Start timer
    time(&time_start);

    // Keep adding to the covariance matrix in n_autotune-sized batches
    for(i_iteration=0;i_iteration<n_autotune;i_iteration++){
      generate_MCMC_chain(MCMC);

      // Build the covariance matrix
      for(i_P=0,k_P=0;i_P<n_P;i_P++){
        P_i_bar_accum[i_P]+=P_chain[i_P];
        for(j_P=0;j_P<n_P;j_P++,k_P++)
          P_ij_bar_accum[k_P]+=P_chain[i_P]*P_chain[j_P];
      }
      n_covariance++;
    }

    // Stop timer
    time(&time_stop);
    time_diff = difftime(time_stop,time_start);

    // Finish the covariance matrix
    for(i_P=0,k_P=0;i_P<n_P;i_P++)
      for(j_P=0;j_P<n_P;j_P++,k_P++)
        V_compute[k_P]=(P_ij_bar_accum[k_P]/(double)n_covariance)-((P_i_bar_accum[i_P]/(double)n_covariance)*(P_i_bar_accum[j_P]/(double)n_covariance));

    // Compute the change in the covariance matrix
    difference=0.;
    for(i_P=0,k_P=0;i_P<n_P;i_P++){
      for(j_P=0;j_P<n_P;j_P++,k_P++){
        difference_i=1e2*(double)fabs(V_compute_last[k_P]-V_compute[k_P])/V_compute_last[k_P];
        difference  =MAX(difference,difference_i);
      }
    }
    memcpy(V_compute_last,V_compute,(size_t)(n_P_squared)*sizeof(double));

    // Report status
    i_autotune++;
    success=100.*(double)(MCMC->n_success)/(double)MCMC->n_propositions;
    SID_log("Iteration #%04d: Success=%5.1lf%% Max change=%5.1f%%",SID_LOG_COMMENT,i_autotune,success,difference);
    time_diff /= (double)n_autotune;
    SID_log("\tMean time for single model call: %.1f", SID_LOG_COMMENT, time_diff);
  }

  // Commit the covariance matrix now that we've computed it
  set_MCMC_covariance(MCMC,V_compute);

  // Clean-up
  SID_free(SID_FARG P_i_bar_accum);
  SID_free(SID_FARG P_ij_bar_accum);
  SID_free(SID_FARG V_compute);
  SID_free(SID_FARG V_compute_last);

  SID_log("Done.",SID_LOG_CLOSE);
}

