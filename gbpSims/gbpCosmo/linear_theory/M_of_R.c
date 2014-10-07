#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double M_of_R(double R,double z,cosmo_info *cosmo){
  double Omega_M;
  double M;
  double rho_bar;
  Omega_M=((double *)ADaPS_fetch((ADaPS *)(cosmo),"Omega_M"))[0];
  rho_bar=Omega_M*rho_crit_z(0.,cosmo);
  M      =FOUR_THIRDS_PI*pow(R,3.)*rho_bar;
  return(M);
}

