#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void fclose_catalog(fp_catalog_info *fp_in){

   sprintf(fp_in->filename_properties_root,"\0");
   sprintf(fp_in->filename_properties_base,"\0");
   sprintf(fp_in->filename_profiles_root,  "\0");
   sprintf(fp_in->filename_profiles_base,  "\0"); 
   if(fp_in->fp_properties!=NULL) fclose(fp_in->fp_properties);
   if(fp_in->fp_profiles!=NULL)   fclose(fp_in->fp_profiles);
   fp_in->fp_properties       =NULL;
   fp_in->fp_profiles         =NULL;
   fp_in->i_file              =0;
   fp_in->n_files             =0;
   fp_in->n_halos_total       =0;
   fp_in->i_halo              =0;
   fp_in->i_halo_start        =0;
   fp_in->i_halo_stop         =0;
   fp_in->n_halos_file        =0;
   fp_in->flag_read_properties=FALSE;
   fp_in->flag_read_profiles  =FALSE;
   fp_in->flag_multifile      =FALSE;

}

