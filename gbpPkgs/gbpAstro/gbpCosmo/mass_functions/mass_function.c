#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_mass_functions.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Compute comoving mass function
double mass_function(double        M_interp,
                     double        z,
                     cosmo_info  **cosmo,
                     int           select_flag){
  int     i;
  int     n_k;
  double  b_z;
  double *lk_P;
  double *ln_Inv_sigma;
  double  sigma_interp;
  interp_info *interp_ln_Inv_sigma;
  double  dlnInvs_dlogM;
  double  Omega_M,rho_o;
  double  derr;
  double  rval;

  // Specify the linear matter power spectrum here
  int mode     =PSPEC_LINEAR_TF;
  int component=PSPEC_ALL_MATTER;

  // Initialize/compute some misc. cosmology things
  sigma_interp=sqrt(power_spectrum_variance(k_of_M(M_interp,0.,*cosmo),z,cosmo,mode,component));
  Omega_M     =((double *)ADaPS_fetch(*cosmo,"Omega_M"))[0];
  rho_o       =Omega_M*rho_crit_z(0.,*cosmo); // Comoving density
  b_z         =linear_growth_factor(z,*cosmo);

  // Compute mass function
  dlnInvs_dlogM=dln_Inv_sigma_dlogM(cosmo,M_interp,z,mode,component);
  rval         =scaled_mass_function(sigma_interp,select_flag)*dlnInvs_dlogM*rho_o/M_interp;
  return(rval);
}

