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

void compute_MCMC_Chi2(MCMC_info *MCMC,double **M,double *Chi2,int *n_DoF){
  int           i_DS;
  int           i_M;
  int           n_M;
  double       *M_target;
  double       *dM_target;
  MCMC_DS_info *current_DS;
  MCMC_DS_info *next_DS;
  (*Chi2)   =0.;
  (*n_DoF)  =0;
  i_DS      =0;
  current_DS=MCMC->DS;
  while(current_DS!=NULL){
    next_DS  =current_DS->next;
    n_M      =current_DS->n_M;
    M_target =current_DS->M_target;
    dM_target=current_DS->dM_target;
    for(i_M=0;i_M<n_M;i_M++){
      (*Chi2)+=pow((M_target[i_M]-M[i_DS][i_M])/dM_target[i_M],2.);
      (*n_DoF)++;
printf("  %le %le %le %le -- %le -- %le\n",current_DS->array[0][i_M],M_target[i_M],dM_target[i_M],M[i_DS][i_M],M[i_DS][i_M],(*Chi2)/(*n_DoF));
    }
    i_DS++;
    current_DS=next_DS;
  }
}

