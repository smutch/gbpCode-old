#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){
  ADaPS     *cosmo;
  double     Omega_M;
  double     Omega_Lambda;
  double     Omega_k;
  double     h_Hubble;
  double     z;
  if(argc<2 || argc>6){
    fprintf(stderr,"\n Syntax: %s z [h_Hubble] [Omega_M] [Omega_Lambda] [Omega_k]\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else
    z=(double)atof(argv[1]);
  init_cosmo_std(&cosmo);
  if(argc>2){
    h_Hubble=(double)atof(argv[2]);
    ADaPS_store(&cosmo,
               "h_Hubble",
               (void *)(&h_Hubble),
               ADaPS_DOUBLE,
               0,NULL);
  }
  if(argc>3){
    Omega_M     =(double)atof(argv[3]);
    Omega_Lambda=1.-Omega_M;
    ADaPS_store(&cosmo,
               "Omega_M",
               (void *)(&Omega_M),
               ADaPS_DOUBLE,
               0,NULL);
    ADaPS_store(&cosmo,
               "Omega_Lambda",
               (void *)(&Omega_Lambda),
               ADaPS_DOUBLE,
               0,NULL);
  }
  if(argc>4){
    Omega_Lambda=(double)atof(argv[4]);
    ADaPS_store(&cosmo,
               "Omega_Lambda",
               (void *)(&Omega_Lambda),
               ADaPS_DOUBLE,
               0,NULL);
  }
  if(argc>5){
    Omega_k=(double)atof(argv[5]);
    ADaPS_store(&cosmo,
               "Omega_k",
               (void *)(&Omega_k),
               ADaPS_DOUBLE,
               0,NULL);
  }
  printf("%10.3lf\n",D_angular(z,cosmo)/M_PER_MPC);
  free_cosmo(&cosmo);
  return(0);
}
