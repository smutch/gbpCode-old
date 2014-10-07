#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

double rho_NFW(double       r,
               double       M_vir,
               double       z,
               int          mode,
               cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double g_c;
  double rho_o;
  double x;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
  rho_o=M_vir*g_c/(FOUR_PI*pow(R_vir/c_vir,3.));
  x    =r*c_vir/R_vir;
  return(rho_o/(x*pow(1.+x,2.)));
}

