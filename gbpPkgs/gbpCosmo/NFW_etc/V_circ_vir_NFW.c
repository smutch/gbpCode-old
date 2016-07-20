#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

double V_circ_vir_NFW(double       M_vir,
                      double       z,
                      int          mode,
                      cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double r;
  double V2_c_characteristic;
  double g_c;
  double r_o;
  double x;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  r_val=V_circ_NFW(R_vir,M_vir,z,mode,cosmo);

  return(r_val);
}

