#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>

double M_gbpCosmo2gbpCosmo(double M,gbpCosmo2gbpCosmo_info *cosmo2cosmo){
  double M_prime=M;

  // Perform scaling if it is given
  if(cosmo2cosmo!=NULL)
     M_prime=M_prime*cosmo2cosmo->s_M;

  return(M_prime);
}

