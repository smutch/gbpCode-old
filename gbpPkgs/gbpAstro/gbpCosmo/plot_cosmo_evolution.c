#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){
  cosmo_info *cosmo;
  double      Omega_M;
  double      Omega_Lambda;
  double      Omega_k;
  double      h_Hubble;
  double      z;
  
  init_cosmo_std(&cosmo);

  for(z=0;z<10.;z+=0.005)
    printf("%10.3le %10.3le %10.3le\n",
           z,
           Omega_z(z,cosmo),
           Delta_vir(z,cosmo));

  free_cosmo(&cosmo);
  return(0);
}
