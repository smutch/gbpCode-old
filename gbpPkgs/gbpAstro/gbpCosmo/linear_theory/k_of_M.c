#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double k_of_M(double M,double z,cosmo_info *cosmo){
  double R;
  double k;
  R=R_of_M(M,z,cosmo);
  k=k_of_R(R);
  return(k);
}

