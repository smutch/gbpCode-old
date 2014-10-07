#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double dln_Inv_sigma_dlogM(cosmo_info **cosmo,
            double       M_interp,
            double       z,
            int          mode,
            int          component){
  double       r_val;
  double      *ln_Inv_sigma;
  double       derr;
  int          i;
  int          n_k;
  size_t  n_k_dim;
  double      *lk_P;
  double      *lM_k;
  interp_info *interp_ln_Inv_sigma;
  char         mode_name[ADaPS_NAME_LENGTH];
  char         component_name[ADaPS_NAME_LENGTH];
  char         interp_name[ADaPS_NAME_LENGTH];

  // Initialize
  pspec_names(mode,component,mode_name,component_name);
  sprintf(interp_name,"ln_Inv_sigma_lnM_%s_%s_interp",mode_name,component_name);
  if(!ADaPS_exist(*cosmo,interp_name))
     init_sigma_M(cosmo,
                  z,
                  mode,
                  component);
  interp_ln_Inv_sigma=(interp_info *)ADaPS_fetch(*cosmo,interp_name);

  double norm =linear_growth_factor(z,*cosmo);
  r_val=take_ln(10.)*
    interpolate_derivative(interp_ln_Inv_sigma,
                 take_ln(M_interp))/norm;
  return(r_val);
}

