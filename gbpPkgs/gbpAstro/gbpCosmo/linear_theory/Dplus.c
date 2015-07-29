#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double Dplus_integrand(double a,void *cosmo){
  return(dDplus_da(a,(cosmo_info *)cosmo));
}

double Dplus(double a,cosmo_info *cosmo){

  // Initialize integral
  int    n_int       =1000;
  double rel_accuracy=1e-8;
  double abs_error;
  double r_val;
  double limit_lo=0.;
  double limit_hi=a;
  gsl_function               integrand;
  gsl_integration_workspace *wspace;
  integrand.function=Dplus_integrand;
  integrand.params  =(void *)cosmo;
  wspace            =gsl_integration_workspace_alloc(n_int);

  // Perform integral
  gsl_integration_qags(&integrand,
                       0.,limit_hi,
                       1e-5,0.,
                       n_int,
                       wspace,
                       &r_val,&abs_error); // use qags for singularity at a=0

  //gsl_integration_qag(&integrand,
  //                     limit_lo,limit_hi,
  //                     0.,rel_accuracy,
  //                     n_int,
  //                     GSL_INTEG_GAUSS61,
  //                     wspace,
  //                     &r_val,&abs_error); // use qags for singularity at a=0

  // Clean-up
  gsl_integration_workspace_free(wspace);

  // Fetch some cosmology stuff
  double Omega_M     =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  double Omega_k     =((double *)ADaPS_fetch(cosmo,"Omega_k"))[0];
  double Omega_Lambda=((double *)ADaPS_fetch(cosmo,"Omega_Lambda"))[0];
  double Ez          =E_z(Omega_M,Omega_k,Omega_Lambda,z_of_a(a));

  // Apply coefficients and return result
  r_val*=2.5*Omega_M*Ez;

  return(r_val);
}

