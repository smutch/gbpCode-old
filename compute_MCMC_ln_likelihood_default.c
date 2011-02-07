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

void compute_MCMC_ln_likelihood_default(MCMC_info *MCMC,double **M,double *P,double *ln_likelihood){
  int           i_DS;
  int           i_M;
  int           n_M;
  double       *M_target;
  double       *dM_target;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;
  (*ln_likelihood)=0.;
  i_DS            =0;
  current_DS      =MCMC->DS;
  while(current_DS!=NULL){
    next_DS  =current_DS->next;
    n_M      =current_DS->n_M;
    M_target =current_DS->M_target;
    dM_target=current_DS->dM_target;
    for(i_M=0;i_M<n_M;i_M++){
      (*ln_likelihood)+=pow((M_target[i_M]-M[i_DS][i_M])/dM_target[i_M],2.);
//fprintf(stderr,"%le %le %le %le\n",M_target[i_M],M[i_DS][i_M],dM_target[i_M],(*ln_likelihood));
    }
    i_DS++;
    current_DS=next_DS;
  }
//if(MCMC->n_map_calls>4) exit(1);
  (*ln_likelihood)*=-0.5;
}

