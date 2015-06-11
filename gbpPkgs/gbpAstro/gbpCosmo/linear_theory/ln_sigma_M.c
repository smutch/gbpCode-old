#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double ln_sigma_M(cosmo_info *cosmo,
                  double      M_interp,
                  double      z,
                  int         mode,
                  int         component){
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         sigma_lnM_name[ADaPS_NAME_LENGTH];
  interp_info *interp_sigma_lnM;
  double       r_val;

  // Initialize
  pspec_names(mode,component,mode_name,component_name);
  sprintf(sigma_lnM_name,"sigma_lnM_%s_%s_interp",mode_name,component_name);
  if(!ADaPS_exist(cosmo,sigma_lnM_name))
     init_sigma_M(&cosmo,
                  mode,
                  component);
  interp_sigma_lnM=(interp_info *)ADaPS_fetch(cosmo,sigma_lnM_name);

  // Perform interpolation
  double norm;
  norm =linear_growth_factor(z,cosmo);
  r_val=take_ln(norm*interpolate(interp_sigma_lnM,take_ln(M_interp)));
  return(r_val);
}

