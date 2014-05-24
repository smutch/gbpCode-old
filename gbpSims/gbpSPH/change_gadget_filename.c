#include <stdio.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <string.h>

void change_gadget_filename(char *filename_root_in,char *filename_root,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename){

  // Determine/set the filename root and path
  char filename_path[MAX_FILENAME_LENGTH];
  strcpy(filename_path,filename_root_in);
  strip_file_root(filename_path);

  if(flag_file_type==0)
     sprintf(filename,"%s/%s_%03d/%s_%03d",filename_path,filename_root,snapshot_number,filename_root,snapshot_number);
  else if(flag_file_type==1)
     sprintf(filename,"%s/%s_%03d",filename_path,filename_root,snapshot_number);
  else if(flag_file_type==2)
     sprintf(filename,"%s/%s",filename_path,filename_root);
  else if(flag_file_type==3)
     sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
  else if(flag_file_type==4)
     sprintf(filename,"%s",filename_root_in);

  if(flag_multifile)
    sprintf(filename,"%s.%d",filename,multifile_number);
}

