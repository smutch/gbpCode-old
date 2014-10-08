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
  double z_lo;
  double z_hi;
  if(argc<3 || argc>4){
    fprintf(stderr,"\n Syntax: %s z [gbpCosmo_file.txt]\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    z_lo=(double)atof(argv[1]);
    z_hi=(double)atof(argv[2]);
  }

  SID_log("Computing angular diameter distance between redshifts z=%lf and z=%lf...",SID_LOG_OPEN,z_lo,z_hi);

  // Initialize cosmology
  ADaPS *cosmo;
  if(argc==2)
     init_cosmo_default(&cosmo);
  else if(argc==3)
     read_gbpCosmo_file(&cosmo,argv[3]);

  SID_log("result=%10.3lf",SID_LOG_COMMENT,D_angular_1to2(z_lo,z_hi,cosmo)/M_PER_MPC);

  // Clean-up
  free_cosmo(&cosmo);

  SID_exit(ERROR_NONE);
}
