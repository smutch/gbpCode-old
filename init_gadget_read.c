#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int init_gadget_read(char *filename_root_in,int snapshot_number,int *flag_multifile,int *flag_file_type,gadget_header_info *header){
  char  filename[MAX_FILENAME_LENGTH];
  int   i_file;
  int   flag_filefound=FALSE;
  int   record_length_in;
  FILE *fp;

  // Determine file format
  for(i_file=0;i_file<5 && !flag_filefound;i_file++){  
    if(i_file==0)
      sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
    else if(i_file==1)
      sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
    else if(i_file==2)
      sprintf(filename,"%s",filename_root_in);
    else if(i_file==3)
      sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
    else if(i_file==4)
      sprintf(filename,"%s",filename_root_in);
    fp=fopen(filename,"r");
    if(fp!=NULL){
      flag_filefound   =TRUE;
      (*flag_multifile)=FALSE;
      (*flag_file_type)=i_file;
    }
    // ... if that doesn't work, check for multi-file
    else{
      strcat(filename,".0");
      fp=fopen(filename,"r");
      if(fp!=NULL){
        flag_filefound   =TRUE;
        (*flag_multifile)=TRUE;
        (*flag_file_type)=i_file;
      }
    }
  }

  if(flag_filefound){
    fread(&record_length_in,sizeof(int),1,fp);
    fread(header,sizeof(gadget_header_info),1,fp);
    fclose(fp);
  }
  return(flag_filefound);
}

