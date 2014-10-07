#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double Dplus_integrand(double  a,
             void   *cosmo){
  return(dDplus_da(a,(cosmo_info *)cosmo));
}
double Dplus(double a,cosmo_info *cosmo){
  int    n_int=1000;
  double h_Hubble;
  double Omega_M;
  double Ez;
  double rel_accuracy=1e-5;
  double abs_error;
  double limit_lo;
  double limit_hi;
  double integral;
  double r_val;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  // Initialize integral
  limit_lo          =0.;
  limit_hi          =a;
  integrand.function=Dplus_integrand;
  integrand.params  =(void *)cosmo;
  wspace            =gsl_integration_workspace_alloc(n_int);

  gsl_integration_qags(&integrand,
             limit_lo,limit_hi,
             0,rel_accuracy,
             n_int,
             wspace,
             &r_val,&abs_error); // use qags for singularity at a=0

  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  Omega_M =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  Ez      =H_z(z_of_a(a),cosmo)/(100.*h_Hubble);

  r_val  *=2.5*Omega_M*Ez;

  // Clean-up
  gsl_integration_workspace_free(wspace);

  return(r_val);
}

