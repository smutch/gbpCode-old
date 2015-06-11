#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double M_of_k(double k,cosmo_info *cosmo){
  double R3     =pow(R_of_k(k),3.);
  double Omega_M=((double *)ADaPS_fetch((ADaPS *)(cosmo),"Omega_M"))[0];
  double rho_bar=Omega_M*rho_crit_z(0,cosmo);
  double M      =FOUR_THIRDS_PI*R3*rho_bar;
  return(M);
}

