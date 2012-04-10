#include <stdio.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void set_gadget_filename(char *filename_root_in,int snapshot_number,int multifile_number,int flag_multifile,int flag_file_type,char *filename){

  if(flag_file_type==0)
    sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
  else if(flag_file_type==1)
    sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
  else if(flag_file_type==2)
    sprintf(filename,"%s",filename_root_in);
  else if(flag_file_type==3)
    sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
  else if(flag_file_type==4)
    sprintf(filename,"%s",filename_root_in);
  if(flag_multifile)
    sprintf(filename,"%s.%d",filename,multifile_number);
}

