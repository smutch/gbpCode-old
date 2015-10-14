#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double t_Hubble_z(double redshift,cosmo_info *cosmo){
  return(1./H_convert(H_z(redshift,cosmo)));
}

