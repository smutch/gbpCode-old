#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void open_gadget_file(char      *filename_root_in,
                      int        snapshot_number,
                      fp_gadget *fp){
  int  i_file;
  char filename[MAX_FILENAME_LENGTH];

  // Determine file format
  for(i_file=0;i_file<3 && !flag_filefound;i_file++){
    if(i_file==0)
      sprintf(filename,"%s/snapshot_%03d/snapshot_%03d",filename_root_in,snapshot_number,snapshot_number);
    else if(i_file==1)
      sprintf(filename,"%s/snapshot_%03d",filename_root_in,snapshot_number);
    else if(i_file==2)
      sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
    fp->fp=fopen(filename,"r");
    if(fp->fp!=NULL){
      fp->flag_filefound=TRUE;
      fp->flag_multifile=FALSE;
      fp->flag_file_type=i_file;
    }
    // ... if that doesn't work, check for multi-file
    else{
      strcat(filename,".0");
      fp->fp=fopen(filename,"r");
      if(fp->fp!=NULL){
        fp->flag_filefound=TRUE;
        fp->flag_multifile=TRUE;
        fp->flag_file_type=i_file;
      }
    }
  }
}

