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

// Compute comoving mass function
double mass_function(double        M_interp,
                     double        z,
                     cosmo_info  **cosmo,
                     int           mode,...){
  va_list vargs;
  va_start(vargs,mode);

  // Initialize/compute some misc. cosmology things
  double sigma_interp=sqrt(power_spectrum_variance(k_of_M(M_interp,*cosmo),z,cosmo,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER));
  double Omega_M     =((double *)ADaPS_fetch(*cosmo,"Omega_M"))[0];
  double rho_o       =Omega_M*rho_crit_z(0,*cosmo); 

  // Compute mass function (eqn 41 of Lukic et al, 2007, for example)
  double dlnInvs_dlogM=dln_Inv_sigma_dlogM(cosmo,M_interp,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);
  double rval         =rho_o*scaled_mass_function(sigma_interp,mode,vargs)*dlnInvs_dlogM/M_interp;
  va_end(vargs);
  return(rval);
}

