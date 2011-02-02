#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void autotune_MCMC_covariance(MCMC_info *MCMC,RNG_info *RNG){
  gsl_matrix    *m;
  gsl_vector    *b;
  int            n_P;
  int            n_P_squared;
  double        *P_new;
  double        *P_last;
  double       **M_new;
  int            n_DS;
  double       **M_last;
  MCMC_DS_info  *current_DS;
  MCMC_DS_info  *next_DS;
  int            i_DS;
  double         threshold;
  double         success;
  int            n_autotune;
  int            n_covariance;
  int            n_success;
  int            i_autotune;
  int            i_iteration;
  double         ln_likelihood_new;
  double         ln_likelihood_last;
  double         ln_Pr_new;
  double         ln_Pr_last;
  int            flag_success;
  double         covariance_threshold;
  double        *P_i_bar_accum;
  double        *P_ij_bar_accum;
  double        *V_compute;
  double        *V_compute_last;
  double         difference;
  double         difference_i;
  int            i_P,j_P,k_P;

  SID_log("Autotuning the covariance matrix (threshold=%6.2lf%%)...",SID_LOG_OPEN|SID_LOG_TIMER,MCMC->covariance_threshold);

  // Initialize parameter arrays
  SID_log("Initizlizing...",SID_LOG_OPEN);
  n_P   =MCMC->n_P;
  n_DS  =MCMC->n_DS;
  P_new =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_last=(double  *)SID_malloc(sizeof(double)  *n_P);
  M_new =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_last=(double  **)SID_malloc(sizeof(double  *)*n_DS);
  current_DS=MCMC->DS;
  i_DS=0;
  while(current_DS!=NULL){
    next_DS     =current_DS->next;
    M_new[i_DS] =(double *)SID_malloc(sizeof(double)*current_DS->n_M);
    M_last[i_DS]=(double *)SID_malloc(sizeof(double)*current_DS->n_M);
    current_DS  =next_DS;
    i_DS++;
  }
  m=MCMC->m;
  b=gsl_vector_calloc(n_P);

  // Initialize covariance matrix
  n_P_squared   =n_P*n_P;
  P_i_bar_accum =(double *)malloc(sizeof(double)*n_P);
  P_ij_bar_accum=(double *)malloc(sizeof(double)*n_P_squared);
  V_compute     =(double *)malloc(sizeof(double)*n_P_squared);
  V_compute_last=(double *)malloc(sizeof(double)*n_P_squared);
  for(i_P=0,k_P=0;i_P<n_P;i_P++){
    P_i_bar_accum[i_P]=0.;
    for(j_P=0;j_P<n_P;j_P++,k_P++){
      P_ij_bar_accum[k_P]=0.;
      V_compute[k_P]     =0.;
    }
  }
  memcpy(V_compute_last,MCMC->V,sizeof(double)*n_P_squared);

  // Initialize everything else
  covariance_threshold=MCMC->covariance_threshold;
  n_autotune          =MCMC->n_autotune_covariance;
  i_autotune          =0;
  difference          =2.*covariance_threshold;
  n_success           =0;
  n_covariance        =0;

  SID_log("Done.",SID_LOG_CLOSE);

  // Integrate the covariance matrix until it convergence
  while(difference>covariance_threshold){

    // Generate the first link to start things off
    if(i_autotune==0){
      memcpy(P_last,MCMC->P_init,(size_t)n_P*sizeof(double));
      while(MCMC->map_P_to_M(P_last,MCMC,M_last)){
        MCMC->first_map_call=FALSE;
        generate_new_MCMC_link(MCMC,P_last,n_P,RNG,m,b,P_new);
        memcpy(P_last,P_new,(size_t)n_P*sizeof(double));
      }
      MCMC->compute_MCMC_ln_likelihood(MCMC,M_last,P_last,&ln_likelihood_last);
    }

    // Keep adding to the covariance matrix in batches
    for(i_iteration=0;i_iteration<n_autotune;i_iteration++){
      generate_new_MCMC_link(MCMC,P_last,n_P,RNG,m,b,P_new);
      while(MCMC->map_P_to_M(P_new,MCMC,M_new)){
        MCMC->first_map_call=FALSE;
        generate_new_MCMC_link(MCMC,P_last,n_P,RNG,m,b,P_new);
      }
      MCMC->compute_MCMC_ln_likelihood(MCMC,M_new,P_new,&ln_likelihood_new);
      ln_Pr_new=MIN(0.0,ln_likelihood_new-ln_likelihood_last);
      if((double)random_number(RNG)<=exp(ln_Pr_new))
        flag_success=TRUE;
      else
        flag_success=FALSE;
      if(flag_success){
        memcpy(P_last,P_new,(size_t)n_P*sizeof(double));
        ln_likelihood_last=ln_likelihood_new;
        ln_Pr_last        =ln_Pr_new;
        n_success++;
      }
      // Build the covariance matrix
      for(i_P=0,k_P=0;i_P<n_P;i_P++){
        P_i_bar_accum[i_P]+=P_last[i_P];
        for(j_P=0;j_P<n_P;j_P++,k_P++)
          P_ij_bar_accum[k_P]+=P_last[i_P]*P_last[j_P];
      }
      n_covariance++;
    }

    // Finish the covariance matrix
    for(i_P=0,k_P=0;i_P<n_P;i_P++)
      for(j_P=0;j_P<n_P;j_P++,k_P++)
        V_compute[k_P]=(P_ij_bar_accum[k_P]/(double)n_covariance)-(P_i_bar_accum[i_P]*P_i_bar_accum[j_P])/((double)n_covariance*(double)n_covariance);

    // Determine the largest change in the covariance matrix and report status
    success=100.*(double)(n_success)/(double)n_covariance;

    // Compute the change in the covariance matrix
    difference=0.;
    for(i_P=0,k_P=0;i_P<n_P;i_P++){
      for(j_P=0;j_P<n_P;j_P++,k_P++){
        difference_i=1e2*(double)fabs(V_compute_last[k_P]-V_compute[k_P])/V_compute_last[k_P];
        difference  =MAX(difference,difference_i);
      }
    }
    SID_log("Iteration #%04d: Success=%5.1lf%% Max change=%5.1f%%",SID_LOG_COMMENT,i_autotune,success,difference);
    memcpy(V_compute_last,V_compute,(size_t)(n_P_squared)*sizeof(double));
    i_autotune++;
  }

  // Commit the covariance matrix now that we've computed it
  set_MCMC_covariance(MCMC,V_compute);

  // Clean-up (remove cmatrix stuff from compute_MCMC too)

  SID_log("Done.",SID_LOG_CLOSE);
}

