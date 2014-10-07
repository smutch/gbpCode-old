#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double sigma2_integrand(double  k,
                        void   *params_in){
  double P_k_tmp;
  double W_k_tmp;
  double r_val;
  sigma2_integrand_params *params;
  params =(sigma2_integrand_params *)params_in;
  P_k_tmp=power_spectrum(k,params->z,&(params->cosmo),params->mode,params->component);
  W_k_tmp=W_k_tophat(k*params->R);
  r_val  =k*k*P_k_tmp*W_k_tmp*W_k_tmp;
  return(r_val);
}

