#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL,NULL);

  // Parse arguments and initialize
  double z;
  if(argc<2 || argc>3){
    fprintf(stderr,"\n Syntax: %s z [gbpCosmo_file.txt]\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else
    z=(double)atof(argv[1]);
  SID_log("Computing cosmology information for z=%.2lf...",SID_LOG_OPEN,z);

  // Initialize cosmology
  ADaPS *cosmo=NULL;
  if(argc==2)
     init_cosmo_default(&cosmo);
  else if(argc==3)
     read_gbpCosmo_file(&cosmo,argv[2]);

  // Output results
  double h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  SID_log("rho_crit     = %13.6le Msol/(Mpc^3)",SID_LOG_COMMENT,rho_crit_z(z,cosmo)*(M_PER_MPC/M_SOL)*M_PER_MPC*M_PER_MPC);
  SID_log("D_angular    = %10.3lf Mpc",         SID_LOG_COMMENT,D_angular(z,cosmo)/M_PER_MPC);
  SID_log("D_luminosity = %10.3lf Mpc",         SID_LOG_COMMENT,D_luminosity(z,cosmo)/M_PER_MPC);
  SID_log("D_comoving   = %10.3lf Mpc",         SID_LOG_COMMENT,D_comove(z,cosmo)/M_PER_MPC);
  SID_log("D_horizon    = %10.3lf Mpc",         SID_LOG_COMMENT,C_VACUUM*deltat_a(&cosmo,0.,a_of_z(z))/M_PER_MPC);
  SID_log("D_V          = %10.3lf Mpc",         SID_LOG_COMMENT,pow((1+z)*D_angular(z,cosmo)*(1+z)*D_angular(z,cosmo)*z*C_VACUUM/H_convert(H_z(z,cosmo)),ONE_THIRD)/M_PER_MPC);
  SID_log("H(z)         = %10.3lf km/s/Mpc",    SID_LOG_COMMENT,H_z       (z, cosmo));
  SID_log("t_age(z)     = %10.3le years",       SID_LOG_COMMENT,t_age_z   (z,&cosmo)/S_PER_YEAR);
  SID_log("t_Hubble(z)  = %10.3le years",       SID_LOG_COMMENT,t_Hubble_z(z, cosmo)/S_PER_YEAR);
  SID_log("t_dyn(z)     = %10.3le years",       SID_LOG_COMMENT,t_dyn_z   (z, cosmo)/S_PER_YEAR);

  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_cosmo(&cosmo);
  SID_exit(ERROR_NONE);
}

