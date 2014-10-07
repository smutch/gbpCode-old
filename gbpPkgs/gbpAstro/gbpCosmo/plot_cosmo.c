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

  // Print syntax
  if(argc<2 || argc>6){
    fprintf(stderr,"\n Syntax: %s z [h_Hubble] [Omega_M] [Omega_Lambda] [Omega_k]\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }

  // Set cosmology
  init_cosmo_std(&cosmo);
  z=(double)atof(argv[1]);
  if(argc>2){
    h_Hubble=(double)atof(argv[2]);
    ADaPS_store(&cosmo,
               (void *)(&h_Hubble),
               "h_Hubble",
               ADaPS_SCALAR_DOUBLE);
  }
  if(argc>3){
    Omega_M     =(double)atof(argv[3]);
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
  if(argc>4){
    Omega_Lambda=(double)atof(argv[4]);
    ADaPS_store(&cosmo,
               (void *)(&Omega_Lambda),
               "Omega_Lambda",
               ADaPS_SCALAR_DOUBLE);
  }
  if(argc>5){
    Omega_k=(double)atof(argv[5]);
    ADaPS_store(&cosmo,
               (void *)(&Omega_k),
               "Omega_k",
               ADaPS_SCALAR_DOUBLE);
  }
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];

  // Output results
  printf("rho_crit     = %13.6le Msol/(Mpc^3)\n",rho_crit_z(z,cosmo)*(M_PER_MPC/M_SOL)*M_PER_MPC*M_PER_MPC);
  printf("D_angular    = %10.3lf Mpc\n",         D_angular(z,cosmo)/M_PER_MPC);
  printf("D_luminosity = %10.3lf Mpc\n",         D_luminosity(z,cosmo)/M_PER_MPC);
  printf("D_comoving   = %10.3lf Mpc\n",         D_comove(z,cosmo)/M_PER_MPC);
  printf("D_horizon    = %10.3lf Mpc\n",         C_VACUUM*deltat_a(&cosmo,0.,a_of_z(z))/M_PER_MPC);
  printf("D_V          = %10.3lf Mpc\n",         pow((1+z)*D_angular(z,cosmo)*(1+z)*D_angular(z,cosmo)*z*C_VACUUM/H_convert(H_z(z,cosmo)),ONE_THIRD)/M_PER_MPC);
  printf("H(z)         = %10.3lf km/s/Mpc\n",    H_z(z,cosmo));
  printf("1/H(z)       = %10.3le years\n",       1./(H_convert(H_z(z,cosmo))*S_PER_YEAR));

  // Clean-up
  free_cosmo(&cosmo);
  SID_exit(ERROR_NONE);
}

