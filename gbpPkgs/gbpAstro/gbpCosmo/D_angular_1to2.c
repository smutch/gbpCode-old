#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){
  ADaPS *cosmo;
  double     Omega_M;
  double     Omega_Lambda;
  double     Omega_k;
  double     h_Hubble;
  double     z_lo,z_hi;
  if(argc<3 || argc>7){
    fprintf(stderr,"\n Syntax: %s z_lo z_hi [h_Hubble] [Omega_M] [Omega_Lambda] [Omega_k]\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    z_lo=(double)atof(argv[1]);
    z_hi=(double)atof(argv[2]);
  }
  init_cosmo_std(&cosmo);
  if(argc>3){
    h_Hubble=(double)atof(argv[3]);
    ADaPS_store(&cosmo,
               (void *)(&h_Hubble),
               "h_Hubble",
               ADaPS_SCALAR_DOUBLE);
  }
  if(argc>4){
    Omega_M     =(double)atof(argv[4]);
    Omega_Lambda=1.-Omega_M;
    ADaPS_store(&cosmo,
               (void *)(&Omega_M),
               "Omega_M",
               ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&cosmo,
               (void *)(&Omega_Lambda),
               "Omega_Lambda",
               ADaPS_SCALAR_DOUBLE);
  }
  if(argc>5){
    Omega_Lambda=(double)atof(argv[5]);
    ADaPS_store(&cosmo,
               (void *)(&Omega_Lambda),
               "Omega_Lambda",
               ADaPS_SCALAR_DOUBLE);
  }
  if(argc>6){
    Omega_k=(double)atof(argv[6]);
    ADaPS_store(&cosmo,
               (void *)(&Omega_k),
               "Omega_k",
               ADaPS_SCALAR_DOUBLE);
  }
  printf("%10.3lf\n",D_angular_1to2(z_lo,z_hi,cosmo)/M_PER_MPC);
  free_cosmo(&cosmo);
  return(0);
}
