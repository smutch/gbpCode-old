#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double Inv_sigma_M(cosmo_info *cosmo,
                   double      M_interp,
                   double      z,
                   int         mode,
                   int         component){
  return(1./(sigma_M(cosmo,M_interp,z,mode,component)));
}

