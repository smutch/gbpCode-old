#include <string.h>
#include <gbpMultifile.h>

int fopen_multifile_nth_file(fp_multifile_info *fp_in,int n){

   // Check that the file number makes sense
   if(fp_in->n_files<n)
      SID_trap_error("Invalid file number (%d) requested for multifile {%s;n_files=%d}.",n,fp_in->filename_base,fp_in->n_files);

   // We can't just jump to the file we want.  We need to keep scaning through them so we know what absolute item range the n'th file represents
   int i_file;
   int r_val=TRUE;

   // Start from the beginning if we are going backwards in the file count
   if(n<fp_in->i_file){
      i_file             =0;
      fp_in->i_item_start=0;
   }
   // We have to scan all files in order, so if we're jumping more than one,
   //    start with the next one first and we'll get there
   else if(n>fp_in->i_file)
      i_file=fp_in->i_file+1;
   else
      i_file=fp_in->i_file;

   // Loop until we get to the file we want (or there's an error)
   while(i_file<=n && r_val){
      char filename_multifile[MAX_FILENAME_LENGTH];

      // Close the file pointers if they are already open
      if(fp_in->fp_multifile!=NULL){
         fclose(fp_in->fp_multifile);
         fp_in->fp_multifile=NULL;
      }

      // Create the filename
      if(fp_in->flag_multifile)
         sprintf(filename_multifile,"%s/%s.%d",fp_in->filename_root,fp_in->filename_base,i_file);
      else{
         if(i_file==0)
            sprintf(filename_multifile,"%s",fp_in->filename_root);
         else
            SID_trap_error("Catalog file identified as non-multi-file {%s} has been accessed as multi-file {requested file=%d}.",ERROR_LOGIC,
                           fp_in->filename_base,i_file);
      }

      // Try to open multifile file
      if(((fp_in->fp_multifile)=fopen(filename_multifile,"r"))==NULL)
         r_val=FALSE;
      else{
         // Read header information
         fread(&(fp_in->i_file),       sizeof(int),1,fp_in->fp_multifile);
         fread(&(fp_in->n_files),      sizeof(int),1,fp_in->fp_multifile);
         fread(&(fp_in->n_items_file), sizeof(int),1,fp_in->fp_multifile);
         fread(&(fp_in->n_items_total),sizeof(int),1,fp_in->fp_multifile);
         // Check that the file number in the file is correct
         if(i_file!=fp_in->i_file)
            SID_trap_error("Invalid file number (ie. %d!=%d) in multifile {%s}.",ERROR_LOGIC,i_file,fp_in->i_file,fp_in->filename_root);
      }

      // Set the absolute start and stop ranges of the item numbers
      if((fp_in->i_file)==0)
         fp_in->i_item_start=0;
      else
         fp_in->i_item_start=fp_in->i_item_stop+1;
      fp_in->i_item_stop=fp_in->i_item_start+fp_in->n_items_file-1;
      fp_in->i_item     =fp_in->i_item_start; // This should always be the item index pointed to by the present state of the file pointers

      i_file++;
   }

   if(!r_val)
      SID_trap_error("Problem encountered opening file {%s/%s;file=%d}",ERROR_LOGIC,fp_in->filename_base,n);

   return(r_val);
}
