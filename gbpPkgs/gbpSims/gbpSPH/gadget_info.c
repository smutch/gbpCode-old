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
  int        snapshot;

  SID_init(&argc,&argv,NULL,NULL);

  // Parse command line
  if(argc!=3){
    fprintf(stderr,"\n syntax: %s filename_root snapshot\n",argv[0]);
    fprintf(stderr," ------\n\n");
    SID_exit(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in, argv[1]);
    snapshot=(int)atoi(argv[2]);
  }

  SID_log("Displaying statistics for GADGET file {%s} ...",SID_LOG_OPEN);

  // Read GADGET file
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_gadget_binary(filename_in,snapshot,&plist,READ_GADGET_DEFAULT);

  // Display info
  display_gadget_info(&plist);

  // Clean-up
  free_plist(&plist);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
