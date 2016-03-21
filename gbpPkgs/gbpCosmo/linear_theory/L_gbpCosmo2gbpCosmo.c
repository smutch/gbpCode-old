#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double L_gbpCosmo2gbpCosmo(double L,gbpCosmo2gbpCosmo_info *cosmo2cosmo){
  double L_prime=L;

  // Perform scaling if it is given
  if(cosmo2cosmo!=NULL)
     L_prime=L_prime*cosmo2cosmo->s_L;

  return(L_prime);
}

