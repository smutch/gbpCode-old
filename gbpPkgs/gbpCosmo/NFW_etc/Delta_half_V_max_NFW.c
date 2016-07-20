#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

// From Alam et al '02
double Delta_half_V_max_NFW(double       M_vir,
                            double       z,
                            int          mode,
                            cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double g_c;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
  r_val=409.*c_vir*c_vir*c_vir*g_c;

  return(r_val);
}

