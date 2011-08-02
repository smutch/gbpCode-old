#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  plist_info plist;
  char       filename_in[256];
  char       filename_out[256];
  int        n_files;

  SID_init(&argc,&argv,NULL);

  // Parse command line
  if(argc!=4){
    fprintf(stderr,"\n syntax: %s filename_in filename_out n_files\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in, argv[1]);
    strcpy(filename_out,argv[2]);
    n_files=atoi(argv[3]);
  }

  SID_log("Rewriting GADGET binary file {%s} in %d parts...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in,n_files);

  // Read GADGET file into data structure
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_gadget_binary(filename_in,&plist,READ_GADGET_MODE_DEFAULT);

  // Rewrite file
  ADaPS_store(&(plist.data),(void *)(&n_files),"n_files",ADaPS_SCALAR_INT);
  write_gadget_binary(filename_out,&plist);

  // Clean-up 
  free_plist(&plist);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}
