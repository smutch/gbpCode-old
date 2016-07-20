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

void compute_MCMC_ln_likelihood_default(MCMC_info *MCMC,double **M,double *P,double *ln_likelihood_DS,int *n_DoF_DS,double *ln_likelihood_all,int *n_DoF_all){
  int           i_DS;
  int           i_M;
  int           n_M;
  double       *M_target;
  double       *dM_target;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;
  (*ln_likelihood_all)=0.;
  (*n_DoF_all)        =0;
  i_DS                =0;
  current_DS          =MCMC->DS;
  while(current_DS!=NULL){
     next_DS               =current_DS->next;
     n_M                   =current_DS->n_M;
     M_target              =current_DS->M_target;
     dM_target             =current_DS->dM_target;
     ln_likelihood_DS[i_DS]=0.;
     n_DoF_DS[i_DS]        =0;
     for(i_M=0;i_M<n_M;i_M++){
        ln_likelihood_DS[i_DS]+=pow((M_target[i_M]-M[i_DS][i_M])/dM_target[i_M],2.);
        n_DoF_DS[i_DS]++;
     }
     (*ln_likelihood_all)  +=ln_likelihood_DS[i_DS];
     ln_likelihood_DS[i_DS]*=-0.5;
     (*n_DoF_all)          +=n_DoF_DS[i_DS];
     n_DoF_DS[i_DS]        -=MCMC->n_P;
     i_DS++;
     current_DS=next_DS;
  }
  (*n_DoF_all)        -=MCMC->n_P;
  (*ln_likelihood_all)*=-0.5;
}

