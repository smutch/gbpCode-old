#include <string.h>
#include <gbpMultifile.h>

void fclose_multifile(fp_multifile_info *fp_in){

   sprintf(fp_in->filename_root,"\0");
   sprintf(fp_in->filename_base,"\0");
   if(fp_in->fp_multifile!=NULL) fclose(fp_in->fp_multifile);
   fp_in->fp_multifile  =NULL;
   fp_in->i_file        =0;
   fp_in->n_files       =0;
   fp_in->n_items_total =0;
   fp_in->i_item        =0;
   fp_in->i_item_start  =0;
   fp_in->i_item_stop   =0;
   fp_in->n_items_file  =0;
   fp_in->flag_multifile=FALSE;
}

