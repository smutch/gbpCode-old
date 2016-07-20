#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double M_sc(double       z,
            cosmo_info **cosmo,
            int          mode,
            int          component){
  int     i;
  int     n_k;
  size_t  n_k_dim;
  double *lM_k;
  double *lk_P;
  double *sigma2;
  interp_info *interp;
  double  delta_sc=1.686;
  double  b_z;
  double  r_val;
  char    mode_name[ADaPS_NAME_LENGTH];
  char    component_name[ADaPS_NAME_LENGTH];
  char    sigma2_name[ADaPS_NAME_LENGTH];
  static double M_sc_last;
  static double z_last=-42.;

  if(z!=z_last){
    // Set/initialize variance
    pspec_names(mode,component,mode_name,component_name);
    sprintf(sigma2_name,"sigma2_k_%s_%s_interp",mode_name,component_name);
    if(!ADaPS_exist(*cosmo,sigma2_name))
      init_power_spectrum_variance(cosmo,mode,component);
    interp   =(interp_info *)ADaPS_fetch(*cosmo,sigma2_name);
    b_z      =linear_growth_factor(z,*cosmo);
    r_val    =bisect_array(interp,delta_sc*delta_sc/(b_z*b_z),1e-4);
    M_sc_last=M_of_k(take_alog10(r_val),*cosmo);
    z_last   =z;
  }

  return(M_sc_last);
}

