#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double rho_crit_z_strip(double redshift,
                        double h_Hubble,
                        double Omega_M,
                        double Omega_Lambda){
  double Hz;
  Hz=H_convert(h_Hubble*1e2*E_z(Omega_M,(1.-Omega_M-Omega_Lambda),Omega_Lambda,redshift));
  return(Hz*Hz/(2.*FOUR_THIRDS_PI*G_NEWTON));
}

