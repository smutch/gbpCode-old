#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double dln_sigma_dlnM(cosmo_info *cosmo,
                      double      M_interp,
                      double      z,
                      int         mode,
                      int         component){

  // Initialize
  char mode_name[ADaPS_NAME_LENGTH];
  char component_name[ADaPS_NAME_LENGTH];
  pspec_names(mode,component,mode_name,component_name);
  if(!ADaPS_exist(cosmo,"ln_sigma_lnM_%s_%s_interp",mode_name,component_name))
     init_sigma_M(&cosmo,
                  mode,
                  component);
  interp_info *interp_ln_sigma=(interp_info *)ADaPS_fetch(cosmo,"ln_sigma_lnM_%s_%s_interp",mode_name,component_name);

  // Perform interpolation
  double norm =linear_growth_factor(z,cosmo);
  double r_val=norm*interpolate_derivative(interp_ln_sigma,take_ln(M_interp));

  return(r_val);
}

