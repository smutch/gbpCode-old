#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double linear_growth_factor(double z,cosmo_info *cosmo){
  // Repeat the integral only if we change 
  //    redshift from the last call.
  static double b_z   =  1.;
  static double z_last=-42.;
  if(z!=z_last){
    double Dplus_a=Dplus(a_of_z(z),cosmo);
    double Dplus_1=Dplus(1.,       cosmo);
    b_z    =Dplus_a/Dplus_1;
    z_last =z;
  }
  return(b_z);
}

