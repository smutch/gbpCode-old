#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double linear_growth_factor(double       redshift,
                            cosmo_info  *cosmo){
  double        Dplus_a;
  double        Dplus_1;
  static double b_z;
  static double z_last=-42.;
  if(redshift!=z_last){
    Dplus_a=Dplus(a_of_z(redshift),cosmo);
    Dplus_1=Dplus(1.,              cosmo);
    b_z    =Dplus_a/Dplus_1;
    z_last =redshift;
  }
  return(b_z);
}

