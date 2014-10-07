#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double E_z(double Omega_M, double Omega_k, double Omega_Lambda, double z){
  double one_plus_z;
  double one_plus_z_sq;
  double one_plus_z_cu;
  double result;
  one_plus_z   =1.+z;
  one_plus_z_sq=one_plus_z*one_plus_z;
  one_plus_z_cu=one_plus_z_sq*one_plus_z;
  result=sqrt(Omega_M*one_plus_z_cu+Omega_k*one_plus_z_sq+Omega_Lambda);
  return(result);
}

