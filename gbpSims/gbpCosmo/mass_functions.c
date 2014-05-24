#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Scaled mass functions
double scaled_mass_function(double sigma,int select_flag){
  double delta_k;
  double rval;
  switch(select_flag){
  case MF_JENKINS: // Jenkins et al. (2001)
    delta_k=1.686;
    rval   =0.315*exp(-pow((double)fabs((float)(take_ln(1./sigma)+0.61)),3.8));
    break;
  case MF_PS: // Press-Schechter (1974)
    delta_k=1.686;
    rval   =sqrt(2./PI)*(delta_k/sigma)*exp(-delta_k*delta_k/(2.*sigma*sigma));
    break;
  case MF_ST: // Sheth-Torman (1999)
    delta_k=1.686;
    rval   =0.3222*sqrt(2.*0.707/PI)*(delta_k/sigma)*exp(-0.707*delta_k*delta_k/(2.*sigma*sigma))*
      (1.+pow(sigma*sigma/(0.707*delta_k*delta_k),0.3));
    break;
  default:
    fprintf(stderr,"ERROR in scaled_mass_function: select_flag=%d unknown.\n",select_flag);
    rval=0.;
    break;
  }
  return rval;
}

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

