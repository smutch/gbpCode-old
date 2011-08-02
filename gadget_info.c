#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  plist_info plist;
  char       filename_in[256];

  SID_init(&argc,&argv,NULL);

  // Parse command line
  if(argc!=2){
    fprintf(stderr,"\n syntax: %s filename_in\n",argv[0]);
    fprintf(stderr," ------\n\n");
    SID_exit(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in, argv[1]);
  }

  SID_log("Displaying statistics for GADGET file {%s} ...",SID_LOG_OPEN);

  // Read GADGET file
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_gadget_binary(filename_in,&plist,READ_GADGET_MODE_DEFAULT);

  // Display info
  display_gadget_info(&plist);

  // Clean-up
  free_plist(&plist);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
