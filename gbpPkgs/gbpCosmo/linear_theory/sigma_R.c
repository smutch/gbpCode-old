#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double sigma_R(cosmo_info **cosmo,
               double       R,
               double       z,
               int          mode,
               int          component){
  return(sigma_M(cosmo,M_of_R(R,*cosmo),z,mode,component));
}

