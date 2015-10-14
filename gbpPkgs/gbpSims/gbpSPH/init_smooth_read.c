#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <string.h>

int init_smooth_read(char *filename_root_in,int snapshot_number,int *flag_multifile,int *flag_file_type,smooth_header_info *header){
  char  filename[MAX_FILENAME_LENGTH];
  char  filename_root[MAX_FILENAME_LENGTH];
  char  filename_path[MAX_FILENAME_LENGTH];
  int   i_file;
  int   flag_filefound=FALSE;
  int   record_length_in;
  FILE *fp;

  // Determine/set the filename root
  strcpy(filename_root,filename_root_in);
  strip_path(filename_root);
  if(strlen(filename_root)==0)
     sprintf(filename_root,"smooth");
  strcpy(filename_path,filename_root_in);
  strip_file_root(filename_path);

  // Determine file format
  for(i_file=0;i_file<5 && !flag_filefound;i_file++){  
     if(i_file==0)
        sprintf(filename,"%s/%s_%03d/%s_%03d",filename_path,filename_root,snapshot_number,filename_root,snapshot_number);
     else if(i_file==1)
        sprintf(filename,"%s/%s_%03d",filename_path,filename_root,snapshot_number);
     else if(i_file==2)
        sprintf(filename,"%s/%s",filename_path,filename_root);
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
     fread(&(header->n_particles_file), sizeof(int),      1,fp);
     fread(&(header->offset),           sizeof(int),      1,fp);
     fread(&(header->n_particles_total),sizeof(long long),1,fp);
     fread(&(header->n_files),          sizeof(int),      1,fp);
     fclose(fp);
  }
  return(flag_filefound);
}

