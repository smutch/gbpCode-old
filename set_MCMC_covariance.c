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

void set_MCMC_covariance(MCMC_info *MCMC,double *V){
  int i_P,j_P;
  if(V==NULL)
    SID_log("Initializing the covariance matrix...",SID_LOG_OPEN);
  else
    SID_log("Updating the covariance matrix...",SID_LOG_OPEN);
  if(MCMC->V==NULL)
    MCMC->V=(double *)SID_malloc(sizeof(double)*MCMC->n_P*MCMC->n_P);
  if(MCMC->m==NULL)
    MCMC->m=gsl_matrix_calloc(MCMC->n_P,MCMC->n_P);
  if(MCMC->b==NULL)
    MCMC->b=gsl_vector_calloc(MCMC->n_P);
  if(V==NULL){
    for(i_P=0;i_P<MCMC->n_P;i_P++){
      for(j_P=0;j_P<MCMC->n_P;j_P++){
        if(i_P==j_P){
          // Matrix decomposition fails if we initialize an element to zero here.
          //   Check for this and start with an different value if true.
          if(MCMC->P_init[i_P]==0.)
            MCMC->V[i_P*MCMC->n_P+j_P]=1e-3*MCMC->P_limit_max[i_P];
          else
            MCMC->V[i_P*MCMC->n_P+j_P]=pow(MCMC->P_init[i_P],2.);
        }
        else
          MCMC->V[i_P*MCMC->n_P+j_P]=0.;
      }
    }
  }
  else
    memcpy(MCMC->V,V,(size_t)(MCMC->n_P*MCMC->n_P)*sizeof(double));
  for(i_P=0;i_P<MCMC->n_P;i_P++)
    for(j_P=0;j_P<MCMC->n_P;j_P++)
      gsl_matrix_set(MCMC->m,i_P,j_P,MCMC->V[i_P*MCMC->n_P+j_P]);
  gsl_linalg_cholesky_decomp(MCMC->m);

  // Set the upper-diagonal to zero to get L where covariance=L*L^T
  for(i_P=0;i_P<MCMC->n_P;i_P++)
    for(j_P=i_P+1;j_P<MCMC->n_P;j_P++)
      gsl_matrix_set(MCMC->m,i_P,j_P,0.);
  SID_log("Done.",SID_LOG_CLOSE);
}

