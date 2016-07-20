#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double R_Delta_z(double      M_Delta,
         double      Delta,
         double      redshift,
         cosmo_info *cosmo){
  double rho_crit;
  rho_crit=rho_crit_z(redshift,cosmo);
  return(pow(M_Delta/(FOUR_THIRDS_PI*Delta*rho_crit),ONE_THIRD));
}

