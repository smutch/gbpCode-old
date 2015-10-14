#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double dDplus_da(double a,cosmo_info *cosmo){
  double Omega_M     =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  double Omega_k     =((double *)ADaPS_fetch(cosmo,"Omega_k"))[0];
  double Omega_Lambda=((double *)ADaPS_fetch(cosmo,"Omega_Lambda"))[0];
  double Ez          =E_z(Omega_M,Omega_k,Omega_Lambda,z_of_a(a));
  return(1.0/pow(a*Ez,3.0));
}

