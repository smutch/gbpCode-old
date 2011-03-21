#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void free_MCMC(MCMC_info *MCMC){
  int           i_P,i_DS;
  int           i_array;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;

  SID_log("Freeing MCMC structure...",SID_LOG_OPEN);

  // Parameter arrays
  for(i_P=0;i_P<MCMC->n_P;i_P++)
    SID_free(SID_FARG MCMC->P_names[i_P]);
  SID_free(SID_FARG MCMC->P_names);
  SID_free(SID_FARG MCMC->problem_name);
  SID_free(SID_FARG MCMC->P_init);
  SID_free(SID_FARG MCMC->P_new);
  SID_free(SID_FARG MCMC->P_last);
  SID_free(SID_FARG MCMC->P_chain);
  SID_free(SID_FARG MCMC->P_limit_min);
  SID_free(SID_FARG MCMC->P_limit_max);
  if(MCMC->n_arrays>0){
    for(i_array=0;i_array<MCMC->n_arrays;i_array++){
      SID_free(SID_FARG MCMC->array[i_array]);
      SID_free(SID_FARG MCMC->array_name[i_array]);
    }
    SID_free(SID_FARG MCMC->array);
    SID_free(SID_FARG MCMC->array_name);
  }

  // Covariance and displacement vector
  if(MCMC->V!=NULL)
    SID_free(SID_FARG MCMC->V);
  if(MCMC->m!=NULL)
    gsl_matrix_free(MCMC->m);
  if(MCMC->b!=NULL)
    gsl_vector_free(MCMC->b);

  // Random number generator
  if(MCMC->RNG!=NULL)
    free_RNG(MCMC->RNG);

  // Dataset arrays
  current_DS=MCMC->DS;
  i_DS=0;
  while(current_DS!=NULL){
    next_DS=current_DS->next;
    SID_free(SID_FARG MCMC->M_new[i_DS]);
    SID_free(SID_FARG MCMC->M_last[i_DS]);
    SID_free(SID_FARG current_DS->M_target);
    SID_free(SID_FARG current_DS->dM_target);
    if(MCMC->n_arrays>0){
      for(i_array=0;i_array<current_DS->n_arrays;i_array++){
        SID_free(SID_FARG current_DS->array[i_array]);
        SID_free(SID_FARG current_DS->array_name[i_array]);
      }
      SID_free(SID_FARG current_DS->array);
      SID_free(SID_FARG current_DS->array_name);
      SID_free(SID_FARG current_DS->D);
    }
    SID_free(SID_FARG current_DS);
    current_DS=next_DS;
    i_DS++;
  }
  SID_free(SID_FARG MCMC->n_M);
  SID_free(SID_FARG MCMC->M_new);
  SID_free(SID_FARG MCMC->M_last);
  SID_free(SID_FARG MCMC->ln_likelihood_DS);
  SID_free(SID_FARG MCMC->ln_likelihood_DS_best);
  SID_free(SID_FARG MCMC->ln_likelihood_DS_peak);
  SID_free(SID_FARG MCMC->n_DoF_DS);
  SID_free(SID_FARG MCMC->n_DoF_DS_best);
  SID_free(SID_FARG MCMC->n_DoF_DS_peak);

  SID_log("Done.",SID_LOG_CLOSE);
}

