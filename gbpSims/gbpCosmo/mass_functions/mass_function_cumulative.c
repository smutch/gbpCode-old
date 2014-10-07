#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_mass_functions.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Compute cumulative comoving mass function
typedef struct mass_function_cumulative_params_struct mass_function_cumulative_params;
struct mass_function_cumulative_params_struct {
  cosmo_info *cosmo; 
  double      z;
  int         mode;
  int         component;
  int         select_flag;
};
double mass_function_cumulative_integrand(double  lM,
                                          void   *params_in);
double mass_function_cumulative_integrand(double  lM,
                                          void   *params_in){
  double M;
  double n;
  double r_val;
  mass_function_cumulative_params *params;
  params=(mass_function_cumulative_params *)params_in;
  M     =take_alog10(lM);
  n     =mass_function(M,params->z,&(params->cosmo),params->mode,params->component,params->select_flag);
  r_val =n;
  return(r_val);
}
double mass_function_cumulative(double       M_interp,
                                double       z,
                                cosmo_info  *cosmo,
                                int          mode,
                                int          component,
                                int          select_flag){
  double  *lk_P;
  int      n_int=5000;
  double   rel_accuracy=1e-3;
  double   limit_lo;
  double   limit_hi;
  double   r_val;
  double   abs_error;
  gsl_integration_workspace       *wspace;
  gsl_function                     integrand;
  mass_function_cumulative_params  params;

  // Initialize integral
  if(!ADaPS_exist(cosmo,"lk_P"))
    init_power_spectrum_TF(&cosmo);
  lk_P              =(double *)ADaPS_fetch(cosmo,"lk_P");
  limit_lo          =take_log10(M_interp);
  limit_hi          =take_log10(M_of_k(take_alog10(lk_P[0]),z,cosmo));
  if(limit_lo<limit_hi){
    integrand.function=mass_function_cumulative_integrand;
    params.z          =z;
    params.cosmo      =cosmo;
    params.mode       =mode;
    params.component  =component;
    params.select_flag=select_flag;
    integrand.params  =(void *)(&params);
    wspace            =gsl_integration_workspace_alloc(n_int);

    // Integrate mass function * halo occupation
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
  else
    r_val=0.;

  return(r_val);

}

