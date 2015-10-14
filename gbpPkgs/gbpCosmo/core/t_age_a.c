#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double t_age_a(double a_in,cosmo_info **cosmo){
  return(deltat_a(cosmo,DELTAT_A_MIN_A,a_in));
}

