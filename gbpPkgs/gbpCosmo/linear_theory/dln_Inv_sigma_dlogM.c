#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double dln_Inv_sigma_dlogM(cosmo_info **cosmo,
                           double       M_interp,
                           int          mode,
                           int          component){

  // Make sure the sigma arrays have been initialized
  char mode_name[ADaPS_NAME_LENGTH];
  char component_name[ADaPS_NAME_LENGTH];
  pspec_names(mode,component,mode_name,component_name);
  if(!ADaPS_exist(*cosmo,"ln_Inv_sigma_lnM_%s_%s_interp",mode_name,component_name))
     init_sigma_M(cosmo,
                  mode,
                  component);

  // Perform interpolation and set return value
  interp_info *interp_ln_Inv_sigma=(interp_info *)ADaPS_fetch(*cosmo,"ln_Inv_sigma_lnM_%s_%s_interp",mode_name,component_name);
  double r_val=take_ln(10.)*(interpolate_derivative(interp_ln_Inv_sigma,take_ln(M_interp)));

  return(r_val);
}

