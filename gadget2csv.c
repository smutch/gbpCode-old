#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  plist_info plist;
  char       filename_in[256];
  char       filename_out[256];

  SID_init(&argc,&argv,NULL);

  /**********************/
  /* Parse command line */
  /**********************/
  if(argc!=2){
    fprintf(stderr,"\n syntax: %s gadget_file\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in, argv[1]);
    strcpy(filename_out,argv[1]);
    strcat(filename_out,".csv");
  }

  SID_log("Converting GADGET file to .csv...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in,filename_out);

  /****************************************/
  /* Read GADGET file into data structure */
  /****************************************/
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  SID_log("Reading GADGET file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);
  read_gadget_binary(filename_in,&plist);
  SID_log("Done.",SID_LOG_CLOSE);

  /********************/
  /* Write ascii file */
  /********************/
  SID_log("Writing .csv file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out);
  write_csv(filename_out,&plist);
  SID_log("Done.",SID_LOG_CLOSE);

  /************/
  /* Clean-up */
  /************/
  free_plist(&plist);

  SID_log("Done.",SID_LOG_CLOSE);

  return(ERROR_NONE);
}
