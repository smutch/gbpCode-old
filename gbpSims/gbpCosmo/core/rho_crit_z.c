#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double rho_crit_z(double redshift, cosmo_info *cosmo){
  double Hz;
  Hz   =H_convert(H_z(redshift,cosmo));
  return(Hz*Hz/(2.*FOUR_THIRDS_PI*G_NEWTON));
}

