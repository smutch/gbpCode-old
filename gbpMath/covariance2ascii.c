#define  _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpMCMC.h>

int main(int argc, char *argv[]){
  char  filename_in1[MAX_FILENAME_LENGTH];
  char  filename_in2[MAX_FILENAME_LENGTH];
  char  filename_in_root[MAX_FILENAME_LENGTH];
  char  filename_out[MAX_FILENAME_LENGTH];
  int   n_iterations_file_total,n_iterations_file_burn;
  double temperature;
  FILE *fp_in1;
  FILE *fp_in2;
  FILE *fp_out;
  double val;
  int    i_P,j_P,n_P,i_C,n_C;

  SID_init(&argc,&argv,NULL,NULL);

  strcpy(filename_in_root,argv[1]);
  sprintf(filename_in2,"%s/results/covariance.dat",  filename_in_root);
  sprintf(filename_in1,"%s/chains/chain_config_%06d.dat",   filename_in_root,0);
  sprintf(filename_out,"%s/results/covariance.ascii",filename_in_root);

  // Skip chain config file's header
  if((fp_in1=fopen(filename_in1,"r"))!=NULL){
    fread(&n_iterations_file_total,sizeof(int),   1,      fp_in1);
    fread(&n_iterations_file_burn, sizeof(int),   1,      fp_in1);
    fread(&temperature,            sizeof(double),1,      fp_in1);
    SID_log("n_iterations_file_total = %d", SID_LOG_COMMENT, n_iterations_file_total);
    SID_log("n_iterations_file_burn  = %d", SID_LOG_COMMENT, n_iterations_file_burn);
    SID_log("temperature             = %.3f", SID_LOG_COMMENT, temperature);
 }
   else
    SID_trap_error("Could not open file {%s}.",ERROR_IO_READ,filename_in1);

  SID_log("Converting covariance matrix file to ascii format...",SID_LOG_OPEN);
  if((fp_in2=fopen(filename_in2,"r"))!=NULL){
    fp_out=fopen(filename_out,"w");
    fread(&n_P,sizeof(int),1,fp_in2);
    SID_log("n_P=%d",SID_LOG_COMMENT,n_P);
    for(i_P=0;i_P<n_P;i_P++){
      for(j_P=0;j_P<n_P;j_P++){
        fread(&val,sizeof(double),1,fp_in1);
        fprintf(fp_out,"%3d %3d %le ",i_P,j_P,val);
        fread(&val,sizeof(double),1,fp_in2);
        fprintf(fp_out,"%le\n",val);
      }
    }
  }
  else
    SID_trap_error("Could not open file {%s}.",ERROR_IO_READ,filename_in2);
  SID_log("Done.",SID_LOG_CLOSE);
  fclose(fp_in1);
  fclose(fp_in2);
  fclose(fp_out);
  SID_exit(ERROR_NONE);
}

