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

void generate_new_MCMC_link(MCMC_info *MCMC,double *P_old,int n_P,double constant,RNG_info *RNG,gsl_matrix *m,gsl_vector *b,double *P_new){
  int         i_P;
  int         j_P;
  int         i_C;
  double      sum;

  // This is used for the very first initialization
  if(constant<=0){
    for(i_P=0;i_P<n_P;i_P++)
      P_new[i_P]=P_old[i_P];
  }
  else{
    // Generate random vector
    for(i_P=0;i_P<n_P;i_P++)
      gsl_vector_set(b,i_P,constant*random_gaussian(RNG));

    // Use a rotated covariance matrix (if it's given)
    if(m!=NULL){
      for(i_P=0;i_P<n_P;i_P++){
        P_new[i_P]=P_old[i_P];
        for(j_P=0;j_P<n_P;j_P++)
          P_new[i_P]+=gsl_matrix_get(m,i_P,j_P)*gsl_vector_get(b,j_P);
      }
    }
    else{
      for(j_P=0;j_P<n_P;j_P++)
        P_new[j_P]=MCMC->P_init[j_P]*(1.+gsl_vector_get(b,j_P));       
    }
    // Enforce parameter limits
    for(i_P=0;i_P<n_P;i_P++){
       if(P_new[i_P]<MCMC->P_limit_min[i_P]) P_new[i_P]=2*(MCMC->P_limit_min[i_P])-P_new[i_P];
       if(P_new[i_P]>MCMC->P_limit_max[i_P]) P_new[i_P]=2*(MCMC->P_limit_max[i_P])-P_new[i_P];
    }
  }
  if(!check_mode_for_flag(MCMC->mode,MCMC_MODE_PARALLEL))
    SID_Bcast(P_new,n_P*sizeof(double),MASTER_RANK);
}

