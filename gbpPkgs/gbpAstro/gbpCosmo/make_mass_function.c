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

  // Initialize cosmology
  cosmo_info *cosmo;
  read_gbpCosmo_file(&cosmo,filename_cosmology);

  SID_log("Constructing mass function between log(M)=%le->%le at z=%lf...",SID_LOG_OPEN,log_M_min,log_M_max,redshift);

  // Decide which parameterization we are going to use
  int select_flag;
  if(!strcmp(paramterization,"JENKINS"))
     select_flag=MF_JENKINS;
  else if(!strcmp(paramterization,"PS"))
     select_flag=MF_PS;
  else if(!strcmp(paramterization,"ST"))
     select_flag=MF_ST;
  else
     SID_trap_error("Invalid parameterization selected {%s}.  Should be {JENKINS,PS or ST}.",ERROR_SYNTAX);

  // Create output filename
  char filename_out[MAX_FILENAME_LENGTH];
  char redshift_text[64];
  float_to_text(redshift,3,redshift_text);
  sprintf(filename_out,"mass_function_z%s.txt",redshift_text);

  // Open file and write header
  FILE *fp_out=NULL;
  fp_out=fopen(filename_out,"w");
  int i_column=1;
  fprintf(fp_out,"# Mass function for %s cosmology at z=%lf\n",filename_cosmology,redshift);
  fprintf(fp_out,"# \n");
  fprintf(fp_out,"# Column (%02d): log M [h^-1 M_sol]\n",i_column++);
  fprintf(fp_out,"#        (%02d): MFctn\n",             i_column++);
  fprintf(fp_out,"#        (%02d): Cumulative MFctn\n",  i_column++);

  // Create the mass function
  for(int i_bin=0;i_bin<n_M_bins;i_bin++){
     double log_M;
     if(i_bin==0)
        log_M=log_M_min;
     else if(i_bin==(n_M_bins-1))
        log_M=log_M_max;
     else
        log_M=log_M_min+(((double)(i_bin))/((double)(n_M_bins-1)))*(log_M_max-log_M_min);
     fprintf(fp_out,"%le %le %le\n",log_M,
                                    mass_function           (take_alog10(log_M),redshift,&cosmo,select_flag),
                                    mass_function_cumulative(take_alog10(log_M),redshift,&cosmo,select_flag));
  }
  fclose(fp_out);
  SID_log("Result written to {%s}.",SID_LOG_COMMENT,filename_out);

  // Clean-up
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

