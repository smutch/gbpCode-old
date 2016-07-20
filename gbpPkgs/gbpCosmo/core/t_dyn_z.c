#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double t_dyn_z(double redshift,cosmo_info *cosmo){
  return(0.1*t_Hubble_z(redshift,cosmo));
}

