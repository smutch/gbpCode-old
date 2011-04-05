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
  free_MCMC_covariance(MCMC);

  // Random number generator
  if(MCMC->RNG!=NULL)
    free_RNG(MCMC->RNG);

  // Dataset arrays
  free_MCMC_arrays(MCMC);
  free_MCMC_DS(MCMC);

  SID_log("Done.",SID_LOG_CLOSE);
}

