#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>

void autotune_MCMC_temperature(MCMC_info *MCMC,RNG_info *RNG){
  gsl_matrix    *m;
  gsl_vector    *b;
  int            n_P;
  double        *P_new;
  double        *P_last;
  double       **M_new;
  int            n_DS;
  double       **M_last;
  MCMC_DS_info  *current_DS;
  MCMC_DS_info  *next_DS;
  int            i_DS;
  double success_target;
  double success_threshold;
  int    n_autotune;
  double success;
  double temperature_hi;
  double temperature_lo;
  double temperature;
  int    flag_raise_range;
  int    i_autotune;
  double ln_likelihood_last;
  int    i_iteration;
  int    n_success;
  double ln_likelihood_new;
  double ln_Pr_new;
  int    flag_success;
  double ln_Pr_last;

  // Initialize parameter arrays
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

  // Initialize everything else
  success_target   =MCMC->success_target;
  success_threshold=MCMC->success_threshold;
  n_autotune       =MCMC->n_autotune_temperature;
  success          =0.;
  temperature_hi   =MCMC->temperature;
  temperature_lo   =0.;
  temperature      =temperature_hi;
  flag_raise_range =TRUE;
  i_autotune       =0;

  // Iterate with a bisection algorythm until convergence
  SID_log("Autotuning the temperature (target=%6.2lf%%, threshold=%6.2lf%%)...",SID_LOG_OPEN|SID_LOG_TIMER,success_target,success_threshold);
  while(fabs(success-success_target)>success_threshold){

    // Perform bisection
    if(i_autotune>0){
      // If the success is too high, the temperature has to go up
      if(success>success_target){
        temperature_lo=temperature;
        if(flag_raise_range)
          temperature_hi*=4.;
      }
      // If the success is too low, the temperature has to go down
      else{
        temperature_hi  =temperature;
        flag_raise_range=FALSE;
      }
      temperature=0.5*(temperature_hi+temperature_lo);
    }
    set_MCMC_temperature(MCMC,temperature);

    // Generate the first link to start things off
    memcpy(P_last,MCMC->P_init,(size_t)n_P*sizeof(double));
    while(MCMC->map_P_to_M(P_last,MCMC,M_last)){
      MCMC->first_map_call=FALSE;
      generate_new_MCMC_link(MCMC,P_last,n_P,RNG,m,b,P_new);
      memcpy(P_last,P_new,(size_t)n_P*sizeof(double));
    }
    MCMC->compute_MCMC_ln_likelihood(MCMC,M_last,P_last,&ln_likelihood_last);

    // Generate a success rate for the current trial temperature
    for(i_iteration=0,n_success=0;i_iteration<n_autotune;i_iteration++){
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
    }
    success=100.*(double)(n_success)/(double)n_autotune;
    i_autotune++;
    SID_log("Iteration #%04d: Temperature=%le Success=%5.1lf%%",SID_LOG_COMMENT,i_autotune,temperature,success);
  }

  // Clean-up
  current_DS=MCMC->DS;
  i_DS=0;
  while(current_DS!=NULL){
    next_DS     =current_DS->next;
    SID_free(SID_FARG M_new[i_DS]);
    SID_free(SID_FARG M_last[i_DS]);
    current_DS  =next_DS;
    i_DS++;
  }
  SID_free(SID_FARG P_new);
  SID_free(SID_FARG P_last);
  SID_free(SID_FARG M_new);
  SID_free(SID_FARG M_last);

  SID_log("Done.",SID_LOG_CLOSE);
}

