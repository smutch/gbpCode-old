#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double sigma2_integrand(double k,void *params_in){
  sigma2_integrand_params *params =(sigma2_integrand_params *)params_in;
  double P_k_tmp  =power_spectrum(k,params->z,&(params->cosmo),params->mode,params->component);
  double k_W_k_tmp=k*W_k_tophat(k*params->R);
  return(P_k_tmp*k_W_k_tmp*k_W_k_tmp);
}

