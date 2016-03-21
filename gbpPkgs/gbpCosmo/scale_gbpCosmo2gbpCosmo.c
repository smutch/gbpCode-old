#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){
   SID_init(&argc,&argv,NULL,NULL);
 
   // Parse arguments and initialize
   char   filename_gbpCosmo_source[MAX_FILENAME_LENGTH];
   char   filename_gbpCosmo_target[MAX_FILENAME_LENGTH];
   double M_min;
   double M_max;
   double z_min;
   double z;
   if(argc!=7){
     fprintf(stderr,"\n Syntax: %s gbpCosmo_file_source.txt gbpCosmo_file_target.txt z_min M_min M_max z\n",argv[0]);
     fprintf(stderr," ------\n\n");
     return(ERROR_SYNTAX);
   }
   else{
     strcpy(filename_gbpCosmo_source,argv[1]);
     strcpy(filename_gbpCosmo_target,argv[2]);
     z_min=(double)atof(argv[3]);
     M_min=(double)atof(argv[4]);
     M_max=(double)atof(argv[5]);
     z    =(double)atof(argv[6]);
   }
   SID_log("Computing scaling of cosmology {%s} to {%s}...",SID_LOG_OPEN,filename_gbpCosmo_source,filename_gbpCosmo_target);
   SID_log("z_min=%lf",             SID_LOG_COMMENT,z_min);
   SID_log("M_min=%le [h^-1 M_sol]",SID_LOG_COMMENT,M_min);
   SID_log("M_max=%le [h^-1 M_sol]",SID_LOG_COMMENT,M_max);
 
   // Initialize cosmology
   ADaPS *cosmo_1=NULL;
   ADaPS *cosmo_2=NULL;
   read_gbpCosmo_file(&cosmo_1,filename_gbpCosmo_source);
   read_gbpCosmo_file(&cosmo_2,filename_gbpCosmo_target);
 
   // Compute scaling
   gbpCosmo2gbpCosmo_info gbpCosmo2gbpCosmo;
   init_gbpCosmo2gbpCosmo(&cosmo_1,&cosmo_2,z_min,M_min*M_SOL,M_max*M_SOL,&gbpCosmo2gbpCosmo);
 
   // Output results
   SID_log("Results:",SID_LOG_OPEN);
   SID_log("s_L=%le",SID_LOG_COMMENT,L_gbpCosmo2gbpCosmo(1.,&gbpCosmo2gbpCosmo));
   SID_log("s_M=%le",SID_LOG_COMMENT,M_gbpCosmo2gbpCosmo(1.,&gbpCosmo2gbpCosmo));
   SID_log("z' =%le",SID_LOG_COMMENT,z_gbpCosmo2gbpCosmo(z, &gbpCosmo2gbpCosmo));
   SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
 
   SID_log("Done.",SID_LOG_CLOSE);
 
   // Clean-up
   free_cosmo(&cosmo_1);
   free_cosmo(&cosmo_2);
   SID_exit(ERROR_NONE);
}

