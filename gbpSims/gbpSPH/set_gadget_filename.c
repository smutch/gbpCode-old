#include <stdio.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <string.h>

void set_gadget_filename(gadget_read_info *fp_gadget,int multifile_number,char *filename){

  // Determine/set the filename root and path
  char filename_root[MAX_FILENAME_LENGTH];
  char filename_path[MAX_FILENAME_LENGTH];
  strcpy(filename_root,fp_gadget->filename_root);
  strip_path(filename_root);
  if(strlen(filename_root)==0)
     sprintf(filename_root,"snapshot");
  strcpy(filename_path,fp_gadget->filename_root);
  strip_file_root(filename_path);

  if(fp_gadget->flag_file_type==0)
     sprintf(filename,"%s/%s_%03d/%s_%03d",filename_path,filename_root,fp_gadget->snapshot_number,filename_root,fp_gadget->snapshot_number);
  else if(fp_gadget->flag_file_type==1)
     sprintf(filename,"%s/%s_%03d",filename_path,filename_root,fp_gadget->snapshot_number);
  else if(fp_gadget->flag_file_type==2)
     sprintf(filename,"%s/%s",filename_path,filename_root);
  else if(fp_gadget->flag_file_type==3)
     sprintf(filename,"%s_%03d",fp_gadget->filename_root,fp_gadget->snapshot_number);
  else if(fp_gadget->flag_file_type==4)
     sprintf(filename,"%s",fp_gadget->filename_root);

  if(fp_gadget->flag_multifile)
    sprintf(filename,"%s.%d",filename,multifile_number);
}

