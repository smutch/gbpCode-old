#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

// Compute Delta_vir(z) following the formula in 
//   Bryan & Norman (ApJ 495, 80, 1998)
double Delta_vir(double      redshift,
                 cosmo_info *cosmo){
  double x;
  double Omega;
  Omega=Omega_z(redshift,cosmo);
  x    =Omega-1.;
  return((18.*PI*PI+82*x-39*x*x)/Omega);
}

