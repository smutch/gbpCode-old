#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

void init_sigma_M(cosmo_info **cosmo,
                  double       z,
                  int          mode,
                  int          component){

  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         d2ln_sigma_name[ADaPS_NAME_LENGTH];
  char         d2ln_Inv_sigma_name[ADaPS_NAME_LENGTH];
  char         sigma_lnM_name[ADaPS_NAME_LENGTH];
  interp_info *interp_ln_sigma;
  interp_info *interp_ln_Inv_sigma;
  interp_info *interp_sigma_lnM;
  pspec_names(mode,component,mode_name,component_name);
  sprintf(d2ln_sigma_name,    "ln_sigma_lnM_%s_%s_interp",    mode_name,component_name);
  sprintf(d2ln_Inv_sigma_name,"ln_Inv_sigma_lnM_%s_%s_interp",mode_name,component_name);
  sprintf(sigma_lnM_name,     "sigma_lnM_%s_%s_interp",       mode_name,component_name);
  if(!ADaPS_exist((*cosmo),d2ln_sigma_name)){
    SID_log("Generating sigma(R) arrays...",SID_LOG_OPEN);
    if(!ADaPS_exist((*cosmo),"lk_P"))
       init_power_spectrum_TF(cosmo);
    int     n_k;
    int     n_k_dim;
    int     i_k;
    double *lk_P;
    double *lM_k;
    double *sigma_lnM;
    double *ln_sigma;
    double *ln_Inv_sigma;
    n_k    =((int   *)ADaPS_fetch((*cosmo),"n_k"))[0];
    lk_P   =(double *)ADaPS_fetch((*cosmo),"lk_P");
    n_k_dim=(size_t)n_k;
    lM_k        =(double *)SID_malloc(sizeof(double)*n_k);
    sigma_lnM   =(double *)SID_malloc(sizeof(double)*n_k);
    ln_sigma    =(double *)SID_malloc(sizeof(double)*n_k);
    ln_Inv_sigma=(double *)SID_malloc(sizeof(double)*n_k);

    SID_log("Computing variances...",SID_LOG_OPEN);
    n_k_dim =(size_t)n_k;
    for(i_k=0;i_k<n_k;i_k++){
      lM_k[i_k]     =take_ln(M_of_k(take_alog10(lk_P[i_k]),z,(*cosmo)));
      sigma_lnM[i_k]=sqrt(power_spectrum_variance(take_alog10(lk_P[i_k]),
                                                  0.,
                                                  cosmo,
                                                  mode,
                                                  component));
      ln_sigma[i_k]    =take_ln(sigma_lnM[i_k]);
      ln_Inv_sigma[i_k]=take_ln(1./sigma_lnM[i_k]);
    }
    SID_log("Done.",SID_LOG_CLOSE);
    SID_log("Initializing interpolation information...",SID_LOG_OPEN);
    init_interpolate(lM_k,
                     ln_sigma,
                     n_k_dim,
                     gsl_interp_cspline,
                     &interp_ln_sigma);
    ADaPS_store_interp(cosmo,
                       (void *)(interp_ln_sigma),
                       d2ln_sigma_name);
    init_interpolate(lM_k,
                     ln_Inv_sigma,
                     (size_t)n_k,
                     gsl_interp_cspline,
                     &interp_ln_Inv_sigma);
    ADaPS_store_interp(cosmo,
                       (void *)(interp_ln_Inv_sigma),
                       d2ln_Inv_sigma_name);
    init_interpolate(lM_k,
                     sigma_lnM,
                     (size_t)n_k,
                     gsl_interp_cspline,
                     &interp_sigma_lnM);
    ADaPS_store_interp(cosmo,
                       (void *)(interp_sigma_lnM),
                       sigma_lnM_name);
    SID_free(SID_FARG sigma_lnM);
    SID_free(SID_FARG ln_sigma);
    SID_free(SID_FARG ln_Inv_sigma);
    SID_log("Done.",SID_LOG_CLOSE);
    SID_log("Done.",SID_LOG_CLOSE);
  }
}


double sigma_M(cosmo_info *cosmo,
               double      M_interp,
               double      z,
               int         mode,
               int         component){
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         sigma_lnM_name[ADaPS_NAME_LENGTH];
  interp_info *interp_sigma_lnM;
  double       r_val;

  // Initialize
  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma_lnM_name,"sigma_lnM_%s_%s_interp",mode_name,component_name);
  if(!ADaPS_exist(cosmo,sigma_lnM_name))
     init_sigma_M(&cosmo,
                  z,
                  mode,
                  component);
  interp_sigma_lnM=(interp_info *)ADaPS_fetch(cosmo,sigma_lnM_name);

  // Perform interpolation
  double norm;
  norm =linear_growth_factor(z,cosmo);
  r_val=norm*interpolate(interp_sigma_lnM,take_ln(M_interp));
  return(r_val);
}

