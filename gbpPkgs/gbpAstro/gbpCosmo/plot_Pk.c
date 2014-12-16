#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  // Parse arguments and initialize
  double z;
  if(argc<2 || argc>3){
    fprintf(stderr,"\n Syntax: %s z [gbpCosmo_file.txt]\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else
    z=(double)atof(argv[1]);

  SID_log("Computing clustering information for z=%.2lf...",SID_LOG_OPEN,z);

  // Initialize cosmology
  cosmo_info *cosmo=NULL;
  if(argc==2)
     init_cosmo_default(&cosmo);
  else if(argc==3)
     read_gbpCosmo_file(&cosmo,argv[2]);

  // Initialize
  int     mode     =PSPEC_LINEAR_TF;
  int     component=PSPEC_ALL_MATTER;
  init_sigma_M(&cosmo,z,mode,component);
  int     n_k     =((int    *)ADaPS_fetch(cosmo,"n_k"))[0];
  double *lk_P    = (double *)ADaPS_fetch(cosmo,"lk_P");
  double  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  SID_log("Done.",SID_LOG_CLOSE);

  // Generate file
  SID_log("Writing table to stdout...",SID_LOG_OPEN);
  double delta_c    =1.686;
  double m_per_mpc_h=M_PER_MPC/h_Hubble;
  printf("# Column (01): k [h Mpc^-1]\n");
  printf("#        (02): P_k [(h^-1 Mpc)^3]\n");
  printf("#        (03): sigma\n");
  for(int i_k=0;i_k<n_k;i_k++){
     double k_P  =take_alog10(lk_P[i_k]);
     double P_k  =power_spectrum(k_P,z,&cosmo,mode,component);
     double sigma=sqrt(power_spectrum_variance(k_P,z,&cosmo,mode,component));
     printf("%10.5le %10.5le %10.5le\n",
            k_P*m_per_mpc_h,
            P_k/pow(m_per_mpc_h,3.),
            sigma);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

