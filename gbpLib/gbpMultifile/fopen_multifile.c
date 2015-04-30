#include <string.h>
#include <gbpMultifile.h>

int fopen_multifile(const char        *filename_root,
                    size_t             data_size,
                    fp_multifile_info *fp_out,
                    ...){
   va_list vargs;
   va_start(vargs,fp_out);
   int   r_val         =TRUE;
   int   flag_filefound=FALSE;

   // Sort out what file format we're working with
   fp_out->fp_multifile=NULL;
   if(SID.I_am_Master){
      int   i_file;
      char  filename_multifile[MAX_FILENAME_LENGTH];

      // Set some filename information
      vsprintf(fp_out->filename_root,filename_root,vargs);
      strcpy(fp_out->filename_base,fp_out->filename_root);
      strip_path(fp_out->filename_base);

      // Try reading a multifile first ...
      sprintf(filename_multifile,"%s/%s.%d",fp_out->filename_root,fp_out->filename_base,0);
      fp_out->fp_multifile=fopen(filename_multifile,"r");
      if(fp_out->fp_multifile==NULL){
         sprintf(filename_multifile,"%s",fp_out->filename_root);
         // ... if we didn't find a multi-file, try reading a single file ...
         fp_out->fp_multifile=fopen(filename_multifile,"r");
         if(fp_out->fp_multifile==NULL){
           r_val=FALSE;
           // Don't report error here so that we can use this routine to check if datasets exist without termination
           //SID_trap_error("Could not open multifile {%s}",ERROR_LOGIC,filename_root);
         }
         // ... we found a single file.  Set flags.
         else
           fp_out->flag_multifile=FALSE;
      }
      // ... we found a multi-file.  Set flags.
      else
         fp_out->flag_multifile=TRUE;

      // Load/set header information
      if(fp_out->fp_multifile!=NULL){
         fread(&(fp_out->i_file),       sizeof(int),1,fp_out->fp_multifile);
         fread(&(fp_out->n_files),      sizeof(int),1,fp_out->fp_multifile);
         fread(&(fp_out->n_items_file), sizeof(int),1,fp_out->fp_multifile);
         fread(&(fp_out->n_items_total),sizeof(int),1,fp_out->fp_multifile);
         fclose(fp_out->fp_multifile);
         fp_out->fp_multifile=NULL;
      }
   }
   SID_Bcast(fp_out,sizeof(fp_multifile_info),MASTER_RANK,SID.COMM_WORLD);

   // Set the data size
   fp_out->data_size=data_size;

   if(r_val){
      // Initialize things by opening the first file.  Even if we want a item
      //   that's deep in the list, we have to scan all the headers (starting
      //   with the first) to find where it is.
      fp_out->fp_multifile=NULL;
      fp_out->i_file      =0;
      fp_out->i_item      =0;
      if(!(r_val=fopen_multifile_nth_file(fp_out,0)))
         SID_trap_error("Error opening multifile file.",ERROR_IO_OPEN);
   }

   va_end(vargs);
   return(r_val);
}

