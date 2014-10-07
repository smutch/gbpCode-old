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
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         d2ln_sigma_name[ADaPS_NAME_LENGTH];
  interp_info *interp_ln_sigma;
  double       r_val;

  // Initialize
  pspec_names(mode,component,mode_name,component_name);
  sprintf(d2ln_sigma_name,"ln_sigma_lnM_%s_%s_interp",mode_name,component_name);
  if(!ADaPS_exist(cosmo,d2ln_sigma_name))
     init_sigma_M(&cosmo,
                  z,
                  mode,
                  component);
  interp_ln_sigma=(interp_info *)ADaPS_fetch(cosmo,d2ln_sigma_name);

  // Perform interpolation
  double norm;
  norm =linear_growth_factor(z,cosmo);
  r_val=norm*
    interpolate_derivative(interp_ln_sigma,
                 take_ln(M_interp));
  return(r_val);
}

