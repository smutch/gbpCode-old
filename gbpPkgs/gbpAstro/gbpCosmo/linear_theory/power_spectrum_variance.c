#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

void init_power_spectrum_variance(cosmo_info **cosmo,double z,int mode,int component){
  int        i;
  int        n_k;
  size_t     n_k_dim;
  double    *lk_P;
  double    *sigma2;
  interp_info *interp;
  double     b_z;
  double     coefficient;
  char       mode_name[ADaPS_NAME_LENGTH];
  char       component_name[ADaPS_NAME_LENGTH];
  char       sigma2_name[ADaPS_NAME_LENGTH];
  char       d2sigma2_name[ADaPS_NAME_LENGTH];
  double     lk_step;
  double     integral;
  int        n_int=5000;
  double     rel_accuracy=5e-5;
  double     limit_lo;
  double     limit_hi;
  double     r_val;
  double     abs_error;
  sigma2_integrand_params    params;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  z=0.;

  // Compute sigma^2 coefficient
  b_z        =1.;
  coefficient=b_z*b_z/(TWO_PI*PI);

  // Make sure k's are initialized
  if(!ADaPS_exist((*cosmo),"n_k")){
    if(mode==PSPEC_LINEAR_TF)
      init_power_spectrum_TF(cosmo);
    else{
      n_k        =200;
      n_k_dim    =(size_t)n_k;
      lk_P       =(double *)SID_malloc(sizeof(double)*n_k_dim);
      lk_P[0]    =take_log10(k_of_R(1e0*M_PER_KPC));
      lk_P[n_k-1]=take_log10(k_of_R(2e4*M_PER_MPC));
      lk_step    =(lk_P[n_k-1]-lk_P[0])/(double)(n_k-1);
      for(i=1;i<n_k-1;i++)
         lk_P[i]=lk_P[i-1]+lk_step;
      ADaPS_store(cosmo,
                  (void *)(&n_k),
                  "n_k",
                  ADaPS_SCALAR_INT);
      ADaPS_store(cosmo,
                  (void *)(lk_P),
                  "lk_P",
                  ADaPS_DEFAULT);
    }
  }
  n_k    =((int   *)ADaPS_fetch((*cosmo),"n_k"))[0];
  lk_P   =(double *)ADaPS_fetch((*cosmo),"lk_P");
  n_k_dim=(size_t)n_k;

  // Initialize integral
  limit_lo          =take_alog10(lk_P[0]);
  limit_hi          =take_alog10(lk_P[n_k-1]);
  integrand.function=sigma2_integrand;
  params.component  =component;
  params.mode       =mode;
  if(mode==PSPEC_LINEAR_TF)
    params.z        =0.;
  else
    params.z        =z;
  params.cosmo      =(*cosmo);
  integrand.params  =(void *)(&params);
  wspace            =gsl_integration_workspace_alloc(n_int);

  // Loop over the whole range of scales
  sigma2=(double *)SID_malloc(sizeof(double)*n_k_dim);
  for(i=0;i<n_k;i++){

    // Perform spherical top-hat integral
    params.R=R_of_k(take_alog10(lk_P[i]));
    gsl_integration_qag(&integrand,
         limit_lo,limit_hi,
         0,rel_accuracy,
         n_int,
         GSL_INTEG_GAUSS61,
         wspace,
         &integral,&abs_error);
    sigma2[i]=coefficient*integral;
  }
  init_interpolate(lk_P,sigma2,(size_t)n_k,gsl_interp_cspline,&interp);


  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma2_name,  "sigma2_k_%s_%s",       mode_name,component_name);
  sprintf(d2sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);
  ADaPS_store(cosmo,
              (void *)(sigma2),
              "sigma2_k_%s_%s",
              ADaPS_DEFAULT,
              mode_name,component_name);
  ADaPS_store_interp(cosmo,
                     (void *)(interp),
                     "sigma2_k_%s_%s_interp",
                     mode_name,component_name);

  // Clean-up
  gsl_integration_workspace_free(wspace);
}

double power_spectrum_variance(double       k_interp,
                               double       redshift,
                               cosmo_info **cosmo,
                               int          mode,
                               int          component){
  int     n_k;
  double *lk_P;
  double *sigma2;
  interp_info *interp;
  double  rval=0.;
  char    mode_name[ADaPS_NAME_LENGTH];
  char    component_name[ADaPS_NAME_LENGTH];
  char    sigma2_name[ADaPS_NAME_LENGTH];
  char    d2sigma2_name[ADaPS_NAME_LENGTH];
  static double z_last=-42.;
  static double norm;

  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma2_name,  "sigma2_k_%s_%s",       mode_name,component_name);
  sprintf(d2sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);

  if(!ADaPS_exist(*cosmo,sigma2_name))
    init_power_spectrum_variance(cosmo,redshift,mode,component);

  n_k     =((int   *)ADaPS_fetch(*cosmo,"n_k"))[0];
  lk_P    =(double *)ADaPS_fetch(*cosmo,"lk_P");
  sigma2  =(double *)ADaPS_fetch(*cosmo,sigma2_name);
  interp  =(interp_info *)ADaPS_fetch(*cosmo,d2sigma2_name);

  norm=pow(linear_growth_factor(redshift,*cosmo),2.);

  rval=norm*interpolate(interp,take_log10(k_interp));

  return(rval);
}

