#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double H_z(double redshift,cosmo_info *cosmo){
  double Omega_M     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_M"))[0];
  double Omega_k     =((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_k"))[0];
  double Omega_Lambda=((double *)ADaPS_fetch((ADaPS *)cosmo,"Omega_Lambda"))[0];
  double h_Hubble    =((double *)ADaPS_fetch((ADaPS *)cosmo,"h_Hubble"))[0];
  return(h_Hubble*1e2*E_z(Omega_M,Omega_k,Omega_Lambda,redshift));
}

