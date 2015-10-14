#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double D_comove_transverse(double z,cosmo_info *cosmo){
  double D_c;
  double D_H;
  double D_c_t;
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  double h_Hubble;
  Omega_M     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_Lambda"))[0];
  h_Hubble    =((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];
  D_c=D_comove(z,cosmo);
  D_H=D_Hubble(h_Hubble);
  if(Omega_k>0.)
    D_c_t=D_H*sinh(sqrt(Omega_k)*D_c/D_H)/sqrt(Omega_k);
  else if(Omega_k==0.)
    D_c_t=D_c;
  else
    D_c_t=D_H*sin(sqrt(fabs(Omega_k))*D_c/D_H)/sqrt(fabs(Omega_k));
  return(D_c_t);
}

