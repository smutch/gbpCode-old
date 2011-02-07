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

void generate_MCMC_parameters(MCMC_info *MCMC){
  int                i_P;
  int                j_P;
  double             factor_i;
  static int         n_P;
  static double      factor;
  static double     *P_new;
  static double     *P_last;
  static double     *P_init;
  static double     *P_limit_min;
  static double     *P_limit_max;
  static gsl_vector *b;
  static gsl_matrix *m;
  static RNG_info   *RNG;

  // Initialize a few things on the first call
  switch(MCMC->first_gen_parameter_call){
    case TRUE:
      n_P        =MCMC->n_P;
      P_new      =MCMC->P_new;
      P_last     =MCMC->P_last;
      P_init     =MCMC->P_init;
      P_limit_min=MCMC->P_limit_min;
      P_limit_max=MCMC->P_limit_max;
      b          =MCMC->b;
      m          =MCMC->m;
      RNG        =MCMC->RNG;
      factor     =2.4/sqrt((double)MCMC->n_M_total);
      MCMC->first_gen_parameter_call=FALSE;
      break;
  }

  // Generate a random displacement vector
  factor_i=MCMC->temperature*factor;
  for(i_P=0;i_P<n_P;i_P++)
    gsl_vector_set(b,i_P,factor_i*random_gaussian(RNG));

  // Use the rotated covariance matrix (if available)
  if(m!=NULL){
    memcpy(P_new,P_last,n_P*sizeof(double));
    for(i_P=0;i_P<n_P;i_P++){
      for(j_P=0;j_P<n_P;j_P++)
        P_new[i_P]+=gsl_matrix_get(m,i_P,j_P)*gsl_vector_get(b,j_P);
    }
  }
  else{
    for(j_P=0;j_P<n_P;j_P++)
      P_new[j_P]=P_init[j_P]*(1.+gsl_vector_get(b,j_P));       
  }

  // Enforce parameter limits
  for(i_P=0;i_P<n_P;i_P++){
     if(P_new[i_P]<P_limit_min[i_P]) P_new[i_P]=2*(P_limit_min[i_P])-P_new[i_P];
     if(P_new[i_P]>P_limit_max[i_P]) P_new[i_P]=2*(P_limit_max[i_P])-P_new[i_P];
  }

  if(!check_mode_for_flag(MCMC->mode,MCMC_MODE_PARALLEL))
    SID_Bcast(P_new,n_P*sizeof(double),MASTER_RANK);
}

