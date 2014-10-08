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
  printf("#        (02): R [h^-1 Mpc] \n");
  printf("#        (03): M [h^-1 M_sol]\n");
  printf("#        (04): V_max [km/s] \n");
  printf("#        (05): P_k [(h^-1 Mpc)^3]\n");
  printf("#        (06): sigma\n");
  printf("#        (07): nu (peak height)\n");
  printf("#        (08): b_BPR\n");
  printf("#        (09): b_TRK\n");
  printf("#        (10): z-space boost Kaiser '87 (applied to b_TRK)\n");
  printf("#        (11): b_TRK total (w/ Kaiser boost)\n");
  printf("#        (12): b_halo_Poole \n");
  printf("#        (13): z-space boost Poole \n");
  printf("#        (14): b_total_Poole \n");
  printf("#        (15): b_halo_Poole        (substructure)\n");
  printf("#        (16): z-space boost Poole (substructure)\n");
  printf("#        (17): b_total_Poole       (substructure)\n");
  for(int i_k=0;i_k<n_k;i_k++){
     double k_P  =take_alog10(lk_P[i_k]);
     double R_P  =R_of_k(k_P);
     double M_R  =M_of_k(k_P,z,cosmo);
     double P_k  =power_spectrum(k_P,0.,&cosmo,mode,component);
     double sigma=sqrt(power_spectrum_variance(k_P,0.,&cosmo,mode,component));
     double V_max=V_max_NFW(M_R,z,NFW_MODE_DEFAULT,&cosmo);
     double nu   =delta_c/sigma;
     double bias =1.;
     if(M_R<1e16*M_SOL){
        printf("%10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le %10.5le\n",
               k_P*m_per_mpc_h,
               R_P/m_per_mpc_h,
               M_R/(M_SOL/h_Hubble),
               V_max*1e-3,
               P_k/pow(m_per_mpc_h,3.),
               sigma,
               nu,
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_BPR),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_TRK),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_TRK|BIAS_MODEL_KAISER_BOOST),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_TRK|BIAS_MODEL_KAISER),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_POOLE_HALO),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_POOLE_ZSPACE),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_POOLE_TOTAL),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_POOLE_SUBSTRUCTURE|BIAS_MODEL_POOLE_HALO),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_POOLE_SUBSTRUCTURE|BIAS_MODEL_POOLE_ZSPACE),
               bias_model(M_R,delta_c,z,&cosmo,BIAS_MODEL_POOLE_SUBSTRUCTURE|BIAS_MODEL_POOLE_TOTAL));
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

