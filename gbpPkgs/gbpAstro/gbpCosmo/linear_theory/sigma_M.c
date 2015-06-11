#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

void init_sigma_M(cosmo_info **cosmo,
                  int          mode,
                  int          component){

  // Initialize array names
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         d2ln_sigma_name[ADaPS_NAME_LENGTH];
  char         d2ln_Inv_sigma_name[ADaPS_NAME_LENGTH];
  char         sigma_lnM_name[ADaPS_NAME_LENGTH];
  pspec_names(mode,component,mode_name,component_name);
  sprintf(d2ln_sigma_name,    "ln_sigma_lnM_%s_%s_interp",    mode_name,component_name);
  sprintf(d2ln_Inv_sigma_name,"ln_Inv_sigma_lnM_%s_%s_interp",mode_name,component_name);
  sprintf(sigma_lnM_name,     "sigma_lnM_%s_%s_interp",       mode_name,component_name);

  // Check if arrays have been initialized yet.  Do so if not.
  if(!ADaPS_exist((*cosmo),d2ln_sigma_name)){
    SID_log("Generating sigma(R) arrays...",SID_LOG_OPEN);
    // Initialize P(k) if it hasn't been done already
    if(!ADaPS_exist((*cosmo),"lk_P"))
       init_power_spectrum_TF(cosmo);

    // Initialize arrays
    int     n_k         =((int   *)ADaPS_fetch((*cosmo),"n_k"))[0];
    double *lk_P        =(double *)ADaPS_fetch((*cosmo),"lk_P");
    double *lM_k        =(double *)SID_malloc(sizeof(double)*n_k);
    double *sigma_lnM   =(double *)SID_malloc(sizeof(double)*n_k);
    double *ln_sigma    =(double *)SID_malloc(sizeof(double)*n_k);
    double *ln_Inv_sigma=(double *)SID_malloc(sizeof(double)*n_k);

    // Populate arrays
    SID_log("Computing z=0 variance arrays...",SID_LOG_OPEN);
    for(int i_k=0;i_k<n_k;i_k++){
      double k_i       =take_alog10(lk_P[i_k]);
      lM_k[i_k]        =take_ln(M_of_k(k_i,(*cosmo)));
      sigma_lnM[i_k]   =sqrt(power_spectrum_variance(k_i,0.,cosmo,mode,component));
      ln_sigma[i_k]    =take_ln(sigma_lnM[i_k]);
      ln_Inv_sigma[i_k]=take_ln(1./sigma_lnM[i_k]);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Initialize and store interpolations
    interp_info *interp_ln_sigma;
    interp_info *interp_ln_Inv_sigma;
    interp_info *interp_sigma_lnM;
    // ln(sigma)(lnM) interpolation
    SID_log("Initializing interpolation information...",SID_LOG_OPEN);
    init_interpolate(lM_k,
                     ln_sigma,
                     (size_t)n_k,
                     gsl_interp_cspline,
                     &interp_ln_sigma);
    ADaPS_store_interp(cosmo,
                       (void *)(interp_ln_sigma),
                       d2ln_sigma_name);
    // ln(1/sigma)(lnM) interpolation
    init_interpolate(lM_k,
                     ln_Inv_sigma,
                     (size_t)n_k,
                     gsl_interp_cspline,
                     &interp_ln_Inv_sigma);
    ADaPS_store_interp(cosmo,
                       (void *)(interp_ln_Inv_sigma),
                       d2ln_Inv_sigma_name);
    // sigma(lnM) interpolation
    init_interpolate(lM_k,
                     sigma_lnM,
                     (size_t)n_k,
                     gsl_interp_cspline,
                     &interp_sigma_lnM);
    ADaPS_store_interp(cosmo,
                       (void *)(interp_sigma_lnM),
                       sigma_lnM_name);
    SID_log("Done.",SID_LOG_CLOSE);

    // Clean-up 
    SID_free(SID_FARG lM_k);
    SID_free(SID_FARG sigma_lnM);
    SID_free(SID_FARG ln_sigma);
    SID_free(SID_FARG ln_Inv_sigma);

    SID_log("Done.",SID_LOG_CLOSE);
  }
}

double sigma_M(cosmo_info *cosmo,
               double      M_interp,
               double      z,
               int         mode,
               int         component){
  // Fetch some needed stuff.  Initialize sigma(M) if needed.
  char mode_name[ADaPS_NAME_LENGTH];
  char component_name[ADaPS_NAME_LENGTH];
  pspec_names(mode,component,mode_name,component_name);
  if(!ADaPS_exist(cosmo,"sigma_lnM_%s_%s_interp",mode_name,component_name))
     init_sigma_M(&cosmo,mode,component);
  interp_info *interp_sigma_lnM=(interp_info *)ADaPS_fetch(cosmo,"sigma_lnM_%s_%s_interp",mode_name,component_name);

  // Perform interpolation
  double norm =linear_growth_factor(z,cosmo);
  double r_val=norm*interpolate(interp_sigma_lnM,take_ln(M_interp));
  return(r_val);
}

