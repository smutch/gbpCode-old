#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  plist_info plist;

  // Set argument defaults
  int     n_out_def=1;
  int     n_x_def  =1;
  int     n_y_def  =1;
  int     n_z_def  =1;
  int     n_x,n_y,n_z,n_out;
  char    filename_in[256],filename_out[256];

  // Initialize
  SID_init(&argc,&argv,NULL);

  // Process user inputs
  if(argc!=7){
    fprintf(stderr,"Syntax: %s n_x n_y n_z n_files_out filename_in filename_out\n",argv[0]);
    fprintf(stderr,"-------\n");
    SID_exit(ERROR_SYNTAX);
  }
  n_x      =(int)atoi(argv[1]);
  n_y      =(int)atoi(argv[2]);
  n_z      =(int)atoi(argv[3]);
  n_out    =(int)atoi(argv[4]);
  strcpy(filename_in, argv[5]);
  strcpy(filename_out,argv[6]);
  
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  // Read GADGET file
  read_gadget_binary(filename_in,&plist,READ_GADGET_MODE_DEFAULT);

  // Store replication info in data structure; needed by write_gadget_binary
  ADaPS_store(&(plist.data),(void *)(&n_out),"n_files",      ADaPS_SCALAR_INT);
  ADaPS_store(&(plist.data),(void *)(&n_x),  "n_x_replicate",ADaPS_SCALAR_INT);
  ADaPS_store(&(plist.data),(void *)(&n_y),  "n_y_replicate",ADaPS_SCALAR_INT);
  ADaPS_store(&(plist.data),(void *)(&n_z),  "n_z_replicate",ADaPS_SCALAR_INT);

  // Rewrite file, replicated the desired number of times
  write_gadget_binary(filename_out,&plist);

  // Clean-up 
  free_plist(&plist);
  SID_exit(ERROR_NONE);
}
