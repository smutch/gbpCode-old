#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double D_luminosity(double z,cosmo_info *cosmo){
  return(D_comove_transverse(z,cosmo)*(1.+z));
}

