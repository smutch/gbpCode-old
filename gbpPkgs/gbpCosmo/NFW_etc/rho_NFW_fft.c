#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>
#include <gbpCosmo_NFW_etc.h>
#include <gsl/gsl_sf_expint.h>

// FFT of NFW profile from White '01
double rho_NFW_fft(double       k,
                   double       M_vir,
                   double       z,
                   int          mode,
                   cosmo_info **cosmo){
  double c_vir;
  double R_vir;
  double r_s;
  double z_k;
  double g_c;
  double rho_o;
  double x1;
  double x2;
  double Ci1;
  double Ci2;
  double Si1;
  double Si2;
  double r_val;
  set_NFW_params(M_vir,z,mode,cosmo,&c_vir,&R_vir);
  r_s  =R_vir/c_vir;
  z_k  =k*r_s;
  g_c  =1./(log(1.+c_vir)-c_vir/(1.+c_vir));
  rho_o=M_vir*g_c/(FOUR_PI*r_s*r_s*r_s);
  x1   =1.+c_vir;
  x2   =c_vir*z_k;
  Ci1=gsl_sf_Ci((1.+c_vir)*z_k);
  Si1=gsl_sf_Si((1.+c_vir)*z_k);
  Ci2=gsl_sf_Ci(           z_k);
  Si2=gsl_sf_Si(           z_k);
  r_val=(FOUR_PI*rho_o*r_s*r_s*r_s/M_vir)*(cos(z_k)*(Ci1-Ci2)+
                       sin(z_k)*(Si1-Si2)-
                       sin(c_vir*z_k)/(z_k*(1.+c_vir)));
  return(r_val);
}

