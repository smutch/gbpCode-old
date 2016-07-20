#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double t_age_z(double redshift,cosmo_info **cosmo){
  return(t_age_a(a_of_z(redshift),cosmo));
}

