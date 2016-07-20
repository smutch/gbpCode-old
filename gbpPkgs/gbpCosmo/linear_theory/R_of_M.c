#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double R_of_M(double M,cosmo_info *cosmo){
  double Omega_M=((double *)ADaPS_fetch((ADaPS *)(cosmo),"Omega_M"))[0];
  double rho_bar=Omega_M*rho_crit_z(0,cosmo);
  double R      =pow(M/(FOUR_THIRDS_PI*rho_bar),ONE_THIRD);
  return(R);
}

