#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

// Cole and Lacey (1996) M(r) profile
double M_r_NFW(double       r,
               double       M_vir,
               double       z,
               int          mode,
               cosmo_info **cosmo){
  double r_o;
  double x;
  double g_c;
  double c_vir;
  double R_vir;
  double r_val=0.;

  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  r_o=R_vir/c_vir;
  x=r/r_o;
  if(x>0.){
    g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
    r_val=M_vir*g_c*(log(1.+x)-x/(1.+x));
  }
  return(r_val);
}

