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

  SID_log("Computing merger rates for z=%.2lf...",SID_LOG_OPEN,z);

  // Initialize cosmology
  cosmo_info *cosmo=NULL;
  if(argc==2)
     init_cosmo_default(&cosmo);
  else if(argc==3)
     read_gbpCosmo_file(&cosmo,argv[2]);

  // Initialize
  int     mode     =PSPEC_LINEAR_TF;
  int     component=PSPEC_ALL_MATTER;
  init_sigma_M(&cosmo,mode,component);
  double  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  SID_log("Done.",SID_LOG_CLOSE);

  // Generate file
  SID_log("Writing table to stdout...",SID_LOG_OPEN);
  double m_per_mpc_h=M_PER_MPC/h_Hubble;
  int    i_column   =1;
  printf("# Column (%02d): zeta\n",    i_column++);
  printf("#        (%02d): Merger rate B/n(zeta) (eqn 10 from Fakhouri & Ma, 2008)\n",i_column++);
  printf("#        (%02d): Merger rate B/n(zeta) (Fit to Millennium data from Fakhouri & Ma, 2010)\n",i_column++);
  double zeta=1.;
  for(int i_zeta=0;i_zeta<20;i_zeta++,zeta*=0.7){
        printf("%10.5le %10.5le %10.5le\n",
               zeta,
               merger_rate_B_n(&cosmo,0,z,1e12*M_SOL,zeta),
               merger_rate_B_n(&cosmo,1,z,1e12*M_SOL,zeta));
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

