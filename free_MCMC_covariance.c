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

void free_MCMC_covariance(MCMC_info *MCMC){
  int           i_DS;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;

  if(MCMC->n_M!=NULL){
    SID_log("Freeing MCMC covariance matrix...",SID_LOG_OPEN);

    SID_free(SID_FARG MCMC->V);
    if(MCMC->m!=NULL){
      gsl_matrix_free(MCMC->m);
      MCMC->m=NULL;
    }
    if(MCMC->b!=NULL){
      gsl_vector_free(MCMC->b);
      MCMC->b=NULL;
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
}


