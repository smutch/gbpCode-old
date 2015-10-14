#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double Ha_Ho(double a, cosmo_info *cosmo){
  double Omega_M;
  double Omega_k;
  double Omega_Lambda;
  Omega_M     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_M"))[0];
  Omega_k     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  Omega_Lambda=((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_Lambda"))[0];
  return(E_z(Omega_M,Omega_k,Omega_Lambda,z_of_a(a))/
         E_z(Omega_M,Omega_k,Omega_Lambda,0.));
}

