#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double dDplus_dz(double z,cosmo_info *cosmo){
  double a=a_of_z(z);
  return(-dDplus_da(a,cosmo)/(a*a)); // = dD/da * da/dz
}

