#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  plist_info plist;
  char       filename_root_in[256];
  char       filename_out[256];
  char       filename[256];
  int        snapshot_number;
  gadget_header_info header;
  int                block_length;
  int  n_files;
  int  i_file;
  int  flag_filefound,flag_multifile,flag_file_type;
  FILE *fp;

  SID_init(&argc,&argv,NULL);

  // Parse command line
  if(argc!=4){
    fprintf(stderr,"\n syntax: %s filename_in_root snapshot_number filename_out\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_root_in,argv[1]);
    snapshot_number=atoi(argv[2]);
    strcpy(filename_out, argv[3]);
  }

  SID_log("Converting GADGET file {%s}, snapshot #%d to an ascii file {%s}...",SID_LOG_OPEN,filename_root_in,snapshot_number,filename_out);
  gadget_read_info fp_gadget;
  if(init_gadget_read(filename_root_in,snapshot_number,&fp_gadget)){
     // Initialize data structure
     init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

     // Read file
     read_gadget_binary(filename_root_in,snapshot_number,&plist,READ_GADGET_DEFAULT);

     // Write ascii file
     write_gadget_ascii(filename_out,&plist);
  }

  // Clean-up
  free_plist(&plist);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
