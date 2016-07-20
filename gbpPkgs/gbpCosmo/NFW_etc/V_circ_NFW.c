#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

// Cole and Lacey (1996) V_c(r) profile
double V_circ_NFW(double       r,
                  double       M_vir,
                  double       z,
                  int          mode,
                  cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double V2_c_characteristic;
  double g_c;
  double r_o;
  double x;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);

  r_o=R_vir/c_vir;
  x=r/r_o;
  if(x>0.){
    V2_c_characteristic=G_NEWTON*M_vir/R_vir;
    g_c=1./(log(1.+c_vir)-c_vir/(1.+c_vir));
    r_val=sqrt(V2_c_characteristic*g_c*(log(1.+x)-x/(1.+x))*c_vir/x);
  }
  return(r_val);
}

