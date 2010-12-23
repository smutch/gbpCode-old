#define  _MAIN
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpMCMC.h>

int main(int argc, char *argv[]){
  char  filename_in[MAX_FILENAME_LENGTH];
  char  filename_in_root[MAX_FILENAME_LENGTH];
  char  filename_out[MAX_FILENAME_LENGTH];
  FILE *fp_in;
  FILE *fp_out;
  double val;
  int    i_P,j_P,n_P,i_C,n_C;

  SID_init(&argc,&argv,NULL);

  strcpy(filename_in_root,argv[1]);
  sprintf(filename_in, "%s/results/covariance.dat",filename_in_root);
  sprintf(filename_out,"%s/results/covariance.ascii",filename_in_root);

  SID_log("Converting covariance matrix file to ascii format...",SID_LOG_OPEN);
  if((fp_in=fopen(filename_in,"r"))!=NULL){
    fp_out=fopen(filename_out,"w");
    fread(&n_P,sizeof(int),1,fp_in);
    SID_log("n_P=%d",SID_LOG_COMMENT,n_P);
    for(i_P=0;i_P<n_P;i_P++){
      for(j_P=0;j_P<n_P;j_P++){
        fread(&val,sizeof(double),1,fp_in);
        fprintf(fp_out,"%3d %3d %le\n",i_P,j_P,val);
      }
    }
  }
  else
    SID_trap_error("Could not open file {%s}.",ERROR_IO_READ,filename_in);
  SID_log("Done.",SID_LOG_CLOSE);
  fclose(fp_in);
  fclose(fp_out);
  SID_exit(ERROR_NONE);
}

