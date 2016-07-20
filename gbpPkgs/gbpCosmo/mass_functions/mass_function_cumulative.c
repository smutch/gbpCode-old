#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_mass_functions.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Compute cumulative comoving mass function
typedef struct mass_function_cumulative_params_struct mass_function_cumulative_params;
struct mass_function_cumulative_params_struct {
  cosmo_info **cosmo; 
  double       z;
  int          mode;
  double      *P;
};
double mass_function_cumulative_integrand(double  lM,
                                          void   *params_in);
double mass_function_cumulative_integrand(double  lM,
                                          void   *params_in){
  mass_function_cumulative_params *params;
  params=(mass_function_cumulative_params *)params_in;
  double M    =take_alog10(lM);
  double r_val=mass_function(M,params->z,params->cosmo,params->mode,params->P);
  return(r_val);
}
double mass_function_cumulative(double       M_interp,
                                double       z,
                                cosmo_info **cosmo,
                                int          mode,...){
  va_list vargs;
  va_start(vargs,mode);

  // Get passed parameters (if mode specifies they are there).
  double *P=NULL;
  if(check_mode_for_flag(mode,MF_PASS_PARAMS))
       P=(double *)va_arg(vargs,double *);

  // Initialize integral
  if(!ADaPS_exist((*cosmo),"lk_P"))
    init_power_spectrum_TF(cosmo);
  double *lk_P    =(double *)ADaPS_fetch((*cosmo),"lk_P");
  double  limit_lo=take_log10(M_interp);
  double  limit_hi=take_log10(M_of_k(take_alog10(lk_P[0]),(*cosmo)));

  // Perform integral
  int    n_int       =5000;
  double rel_accuracy=1e-4;
  double r_val       =0.;
  if(limit_lo<limit_hi){
     double                           abs_error;
     mass_function_cumulative_params  params;
     gsl_function                     integrand;
     gsl_integration_workspace       *wspace;
     params.z    =z;
     params.cosmo=cosmo;
     params.mode =mode;
     params.P    =P;
     integrand.function=mass_function_cumulative_integrand;
     integrand.params  =(void *)(&params);
     wspace            =gsl_integration_workspace_alloc(n_int);

     // Integrate mass function 
     gsl_integration_qag(&integrand,
                         limit_lo,limit_hi,
                         0,rel_accuracy,
                         n_int,
                         GSL_INTEG_GAUSS61,
                         wspace,
                         &r_val,&abs_error);

     // Clean-up
     gsl_integration_workspace_free(wspace);
  }

  va_end(vargs);
  return(r_val);
}

