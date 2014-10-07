#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double D_comove_integrand(double z, void *cosmo){
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double Ez;
  double r_val;
  Omega_M     =((double *)ADaPS_fetch((cosmo_info *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((cosmo_info *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((cosmo_info *)cosmo,"Omega_Lambda"))[0];
  if(z<=0.)
    r_val=1.;
  else{
    Ez   =E_z(Omega_M,Omega_k,Omega_Lambda,z);
    r_val=1.0/Ez;
  }
  return(r_val);
}
double D_comove(double z,cosmo_info *cosmo){
  int    n_int=1000;
  double h_Hubble;
  double rel_accuracy=1e-5;
  double abs_error;
  double limit_lo;
  double limit_hi;
  double integral;
  double r_val;
  gsl_integration_workspace *wspace;
  gsl_function               integrand;

  if (z<=0.)
    return(0.0);

  h_Hubble=((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];

  // Initialize integral
  limit_lo          =0.;
  limit_hi          =z;
  integrand.function=D_comove_integrand;
  integrand.params  =(void *)cosmo;
  wspace            =gsl_integration_workspace_alloc(n_int);

  gsl_integration_qag(&integrand,
            limit_lo,limit_hi,
            0,rel_accuracy,
            n_int,
            GSL_INTEG_GAUSS61,
            wspace,
            &integral,&abs_error);
  r_val=D_Hubble(h_Hubble)*integral;

  // Clean-up
  gsl_integration_workspace_free(wspace);

  return(D_Hubble(h_Hubble)*integral);
}

