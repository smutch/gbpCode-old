#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <string.h>

int init_gadget_read(char *filename_root_in,int snapshot_number,gadget_read_info *fp_gadget){
  // Determine/set the filename root
  int flag_filefound=FALSE;
  if(SID.I_am_Master){
     char  filename[MAX_FILENAME_LENGTH];
     char  filename_root[MAX_FILENAME_LENGTH];
     char  filename_path[MAX_FILENAME_LENGTH];
     int   record_length_in;
     FILE *fp;
     strcpy(filename_root,filename_root_in);
     strip_path(filename_root);
     if(strlen(filename_root)==0)
        sprintf(filename_root,"snapshot");
     strcpy(filename_path,filename_root_in);
     strip_file_root(filename_path);

     // Store inputs
     sprintf(fp_gadget->filename_root,"%s",filename_root_in);
     fp_gadget->snapshot_number=snapshot_number;

     // Determine file format
     for(int i_type=0;i_type<5 && !flag_filefound;i_type++){  
        if(i_type==0)
           sprintf(filename,"%s/%s_%03d/%s_%03d",filename_path,filename_root,snapshot_number,filename_root,snapshot_number);
        else if(i_type==1)
           sprintf(filename,"%s/%s_%03d",filename_path,filename_root,snapshot_number);
        else if(i_type==2)
           sprintf(filename,"%s/%s",filename_path,filename_root);
        else if(i_type==3)
           sprintf(filename,"%s_%03d",filename_root_in,snapshot_number);
        else if(i_type==4)
           sprintf(filename,"%s",filename_root_in);
        fp=fopen(filename,"r");
        if(fp!=NULL){
           flag_filefound             =TRUE;
           (fp_gadget->flag_multifile)=FALSE;
           (fp_gadget->flag_file_type)=i_type;
        }
        // ... if that doesn't work, check for multi-file
        else{
           strcat(filename,".0");
           fp=fopen(filename,"r");
           if(fp!=NULL){
              flag_filefound             =TRUE;
              (fp_gadget->flag_multifile)=TRUE;
              (fp_gadget->flag_file_type)=i_type;
           }
        }
     }

     // Initialize some flags
     fp_gadget->first_select_call=TRUE;
     fp_gadget->first_action_call=TRUE;

     // Read the first header if we've found a file ok
     if(flag_filefound){
        fread_verify(&record_length_in,   sizeof(int),               1,fp);
        fread_verify(&(fp_gadget->header),sizeof(gadget_header_info),1,fp);
        fclose(fp);
     }

     // Count all the particles in the snapshot
     fp_gadget->n_particles=0;
     for(int i_type=0;i_type<N_GADGET_TYPE;i_type++){
        fp_gadget->n_all[i_type] =(fp_gadget->header.n_all_lo_word[i_type])+(((size_t)(fp_gadget->header.n_all_hi_word[i_type]))<<32);
        fp_gadget->n_particles  +=fp_gadget->n_all[i_type];
     }
  }
  SID_Bcast(fp_gadget,      sizeof(gadget_read_info),MASTER_RANK,SID.COMM_WORLD);
  SID_Bcast(&flag_filefound,sizeof(int),             MASTER_RANK,SID.COMM_WORLD);

  return(flag_filefound);
}

