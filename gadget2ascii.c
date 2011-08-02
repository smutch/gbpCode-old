#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
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
    strcpy(filename_out,    argv[3]);
  }

  SID_log("Converting GADGET file {%s}, snapshot #%d to an ascii file {%s}...",SID_LOG_OPEN,filename_root_in,snapshot_number,filename_out);

  // Initialize data structure
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  // Read GADGET header to get the number of files
  // Determine file format
  n_files=1;
  for(i_file=0,flag_filefound=FALSE;i_file<=3 && !flag_filefound;i_file++){  
    if(i_file==0)
      sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
    else if(i_file==1)
      sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
    else if(i_file==2)
      sprintf(filename,"%s",filename_root_in,snapshot_number);
    else if(i_file==3)
      sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
    fp=fopen(filename,"r");
    if(fp!=NULL){
      flag_filefound=TRUE;
      flag_multifile=FALSE;
      flag_file_type=i_file;
      fread(&block_length,sizeof(int),               1,fp);
      fread(&header,      sizeof(gadget_header_info),1,fp);
      n_files=header.n_files;
      fclose(fp);
    }
    // ... if that doesn't work, check for multi-file
    else{
      strcat(filename,".0");
      fp=fopen(filename,"r");
      if(fp!=NULL){
        flag_filefound=TRUE;
        flag_multifile=TRUE;
        flag_file_type=i_file;
        fread(&block_length,sizeof(int),               1,fp);
        fread(&header,      sizeof(gadget_header_info),1,fp);
        n_files=header.n_files;
        fclose(fp);
      }
    }
  }
  if(flag_filefound){
    for(i_file=0;i_file<n_files;i_file++){
      // Construct file name
      if(flag_file_type==0)
        sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
      else if(flag_file_type==1)
        sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
      else if(flag_file_type==2)
        sprintf(filename,"%s",filename_root_in,snapshot_number);
      else if(flag_file_type==3)
        sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
      if(flag_multifile)
        sprintf(filename,"%s.%d",filename,i_file);

      // Read file
      read_gadget_binary(filename,&plist,READ_GADGET_MODE_DEFAULT);

      // Write ascii file
      write_ascii(filename_out,&plist,i_file);
    }
  }
  // Clean-up
  free_plist(&plist);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
