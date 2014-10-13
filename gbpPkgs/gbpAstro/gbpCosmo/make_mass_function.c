#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  char   filename_cosmology[MAX_FILENAME_LENGTH];
  char   paramterization[MAX_FILENAME_LENGTH];
  double log_M_min    =atof(argv[1]);
  double log_M_max    =atof(argv[2]);
  int    n_M_bins     =atoi(argv[3]);
  double redshift     =atof(argv[4]);
  strcpy(filename_cosmology,argv[5]);
  strcpy(paramterization,   argv[6]);

  SID_log("Constructing mass function between log(M)=%5.3lf->%5.3lf at z=%5.3lf...",SID_LOG_OPEN,log_M_min,log_M_max,redshift);

  // Initialize cosmology
  cosmo_info *cosmo=NULL;
  read_gbpCosmo_file(&cosmo,filename_cosmology);

  // Decide which parameterization we are going to use
  int  select_flag;
  char mfn_text[32];
  if(!strcmp(paramterization,"JENKINS")){
     sprintf(mfn_text,"Jenkins");
     select_flag=MF_JENKINS;
  }
  else if(!strcmp(paramterization,"PS")){
     sprintf(mfn_text,"Press & Schechter");
     select_flag=MF_PS;
  }
  else if(!strcmp(paramterization,"ST")){
     sprintf(mfn_text,"Sheth & Tormen");
     select_flag=MF_ST;
  }
  else
     SID_trap_error("Invalid parameterization selected {%s}.  Should be {JENKINS,PS or ST}.",ERROR_SYNTAX,paramterization);

  // Create output filename
  char  filename_out[MAX_FILENAME_LENGTH];
  char  redshift_text[64];
  char *cosmology_name=(char *)ADaPS_fetch(cosmo,"name");
  float_to_text(redshift,3,redshift_text);
  sprintf(filename_out,"mass_function_z%s_%s_%s.txt",redshift_text,cosmology_name,paramterization);

  // Open file and write header
  FILE *fp_out=NULL;
  fp_out=fopen(filename_out,"w");
  int i_column=1;
  fprintf(fp_out,"# Mass function (%s) for %s cosmology at z=%lf\n",mfn_text,filename_cosmology,redshift);
  fprintf(fp_out,"# \n");
  fprintf(fp_out,"# Column (%02d): log M [h^-1 M_sol]\n",                              i_column++);
  fprintf(fp_out,"#        (%02d): Mass function [(h^{-1} Mpc]^{-3} per dlogM]\n",     i_column++);
  fprintf(fp_out,"#        (%02d): Cumulative Mass function(>M) [(h^{-1} Mpc]^{-3}]\n",i_column++);

  // Create the mass function
  SID_log("Writing results to {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out);
  pcounter_info pcounter;
  SID_init_pcounter(&pcounter,n_M_bins,10);
  double h_Hubble   =((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  double mass_factor=M_SOL/h_Hubble;
  double MFct_factor=pow(M_PER_MPC,3.0);
  for(int i_bin=0;i_bin<n_M_bins;i_bin++){
     double log_M;
     if(i_bin==0)
        log_M=log_M_min;
     else if(i_bin==(n_M_bins-1))
        log_M=log_M_max;
     else
        log_M=log_M_min+(((double)(i_bin))/((double)(n_M_bins-1)))*(log_M_max-log_M_min);
     fprintf(fp_out,"%le %le %le\n",log_M,
                                    MFct_factor*mass_function           (mass_factor*take_alog10(log_M),redshift,&cosmo,select_flag),
                                    MFct_factor*mass_function_cumulative(mass_factor*take_alog10(log_M),redshift,&cosmo,select_flag));
     SID_check_pcounter(&pcounter,i_bin);
  }
  fclose(fp_out);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

