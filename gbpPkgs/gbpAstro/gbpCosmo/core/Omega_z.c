#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_core.h>

double Omega_z(double      redshift,
               cosmo_info *cosmo){
  double Omega_M,Omega_k,Omega_Lambda;
  double Ez,one_plus_z_cube;
  Omega_M        =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  Omega_k        =((double *)ADaPS_fetch(cosmo,"Omega_k"))[0];
  Omega_Lambda   =((double *)ADaPS_fetch(cosmo,"Omega_Lambda"))[0];
  Ez             =E_z(Omega_M,Omega_k,Omega_Lambda,redshift);
  one_plus_z_cube=(1.+redshift)*(1.+redshift)*(1.+redshift);
  return(Omega_M*one_plus_z_cube/(Ez*Ez));
}

