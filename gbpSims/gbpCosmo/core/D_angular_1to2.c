#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double D_angular_1to2(double z_1,double z_2,cosmo_info *cosmo){
  double D_c_t_1;
  double D_c_t_2;
  double D_H;
  double h_Hubble;
  double Omega_k;
  h_Hubble=((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];
  Omega_k =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  D_c_t_1 =D_comove_transverse(z_1,cosmo);
  D_c_t_2 =D_comove_transverse(z_2,cosmo);
  D_H     =D_Hubble(h_Hubble);
  return((D_c_t_2*sqrt(1.+Omega_k*pow(D_c_t_1/D_H,2.0))-
     D_c_t_1*sqrt(1.+Omega_k*pow(D_c_t_2/D_H,2.0)))/(1.+z_2));
}

