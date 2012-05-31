#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

int fopen_nth_catalog_file(fp_catalog_info *fp_in,int n){

   // Check that the file number makes sense
   if(fp_in->n_files<(n+1))
      SID_trap_error("Invalid file number (%d) requested for catalog {%s;n_files=%d}.",n,fp_in->filename_properties_base,fp_in->n_files);

   // We can't just jump to the file we want.  We need to keep scaning through them so we know what absolute halo range the n'th file represents
   int i_file=0;
   int r_val=FALSE;
   i_file=fp_in->i_file;

   // Start from the beginning if we are going backwards in the file count
   if(n<i_file){
      fp_in->i_file      =0;
      fp_in->i_halo_start=0;
      fp_in->i_halo_stop =0;
   }

   // Loop until we get to the file we want (or there's an error)
   while(i_file<=n && !r_val){
      char filename_properties[MAX_FILENAME_LENGTH];
      char filename_profiles[MAX_FILENAME_LENGTH];

      // Close the file pointers if they are already open
      if(fp_in->fp_properties!=NULL){
         fclose(fp_in->fp_properties);
         fp_in->fp_properties=NULL;
      }
      if(fp_in->fp_profiles!=NULL){
         fclose(fp_in->fp_profiles);
         fp_in->fp_profiles=NULL;
      }

      // Create the filename
      if(fp_in->flag_multifile){
         sprintf(filename_properties,"%s/%s.%d",fp_in->filename_properties_root,fp_in->filename_properties_base,i_file);
         sprintf(filename_profiles,  "%s/%s.%d",fp_in->filename_profiles_root,  fp_in->filename_profiles_base,  i_file);
      }
      else{
         if(i_file==0){
            sprintf(filename_properties,"%s",fp_in->filename_properties_root);
            sprintf(filename_profiles,  "%s",fp_in->filename_profiles_root);
         }
         else
            SID_trap_error("Catalog file identified as non-multi-file {%s} has been accessed as multi-file {requested file=%d}.",ERROR_LOGIC,
                           fp_in->filename_properties_base,i_file);
      }

      // Try to open properties file
      if(fp_in->flag_read_properties){
         if(((fp_in->fp_properties)=fopen(filename_properties,"r"))==NULL)
           r_val=TRUE;
         else{
            // Read header information
            fread(&(fp_in->i_file),       sizeof(int),1,fp_in->fp_properties);
            fread(&(fp_in->n_files),      sizeof(int),1,fp_in->fp_properties);
            fread(&(fp_in->n_halos_file), sizeof(int),1,fp_in->fp_properties);
            fread(&(fp_in->n_halos_total),sizeof(int),1,fp_in->fp_properties);
            // Check that the file number in the file is correct
            if(i_file!=fp_in->i_file)
               SID_trap_error("Invalid file number (ie. %d!=%d) in catalog {%s}.",ERROR_LOGIC,i_file,fp_in->i_file,fp_in->filename_properties_root);
            i_file++;
         }
      }

      // Try to open profiles file
      if(fp_in->flag_read_profiles){
         if(((fp_in->fp_profiles)=fopen(filename_profiles,"r"))==NULL)
           r_val=TRUE;
         else{
            // Read header information
            fread(&(fp_in->i_file),       sizeof(int),1,fp_in->fp_profiles);
            fread(&(fp_in->n_files),      sizeof(int),1,fp_in->fp_profiles);
            fread(&(fp_in->n_halos_file), sizeof(int),1,fp_in->fp_profiles);
            fread(&(fp_in->n_halos_total),sizeof(int),1,fp_in->fp_profiles);
            // Check that the file number in the file is correct
            if(i_file!=fp_in->i_file)
               SID_trap_error("Invalid file number (ie. %d!=%d) in catalog {%s}.",ERROR_LOGIC,i_file,fp_in->i_file,fp_in->filename_profiles_root);
            i_file++;
         }
      }

      // Set the absolute start and stop ranges of the halo numbers
      if(fp_in->i_file==0){
         fp_in->i_halo_start=0;
         fp_in->i_halo_stop =fp_in->n_halos_file-1;
      }
      else{
         fp_in->i_halo_start=fp_in->i_halo_start+fp_in->n_halos_file;
         fp_in->i_halo_stop =fp_in->i_halo_stop +fp_in->n_halos_file;
      }
      fp_in->i_halo=fp_in->i_halo_start; // This should always be the halo index pointed to by the present state of the file pointers

      i_file++;
   }

   return(r_val);
}

int fopen_catalog(char            *filename_catalog_root,
                  int              snapshot_number,
                  int              mode,
                  fp_catalog_info *fp_out){
   char  filename_root[MAX_FILENAME_LENGTH];
   int   r_val=FALSE;
   int   flag_filefound=FALSE;
   char  group_text_prefix[8];

   // Decide whether we're opening a group or subgroup catalog
   int flag_read_groups   =FALSE;
   int flag_read_subgroups=FALSE;
   if(check_mode_for_flag(mode,READ_CATALOG_GROUPS)){
     flag_read_groups   =TRUE;
     flag_read_subgroups=FALSE;
     sprintf(group_text_prefix,"");
   }
   if(check_mode_for_flag(mode,READ_CATALOG_SUBGROUPS)){
     if(flag_read_groups)
        SID_trap_error("You can't open both groups and subgroups at the same time in fopen_catalog().",ERROR_LOGIC);
     flag_read_groups   =FALSE;
     flag_read_subgroups=TRUE;
     sprintf(group_text_prefix,"sub");
   }
   if(!flag_read_groups && !flag_read_subgroups)
      SID_trap_error("You must specify either READ_CATALOG_GROUPS or READ_CATALOG_SUBGROUPS in mode for fopen_catalog().",ERROR_LOGIC);

   // Decide if we are reading properties
   if(check_mode_for_flag(mode,READ_CATALOG_PROPERTIES))
      fp_out->flag_read_properties=TRUE;
   else
      fp_out->flag_read_properties=FALSE;

   // Decide if we are reading profiles
   if(check_mode_for_flag(mode,READ_CATALOG_PROFILES))
      fp_out->flag_read_profiles=TRUE;
   else
      fp_out->flag_read_profiles=FALSE;

   // Make sure we're reading something
   if(fp_out->flag_read_properties==FALSE && fp_out->flag_read_profiles==FALSE)
      SID_trap_error("Neither properties nor profiles are selected for reading in fopen_catalog().",ERROR_LOGIC);

   // Set snapshot number
   fp_out->snap_num=snapshot_number;

   // Sort out what file format we're working with
   fp_out->fp_properties=NULL;
   fp_out->fp_profiles  =NULL;
   if(SID.I_am_Master){
      int   i_file;
      char  filename_properties[MAX_FILENAME_LENGTH];
      char  filename_profiles[MAX_FILENAME_LENGTH];

      // Set some filename information
      sprintf(fp_out->filename_properties_root,"%s_%03d.catalog_%sgroups_properties",filename_catalog_root,snapshot_number,group_text_prefix);
      sprintf(fp_out->filename_profiles_root,  "%s_%03d.catalog_%sgroups_profiles",  filename_catalog_root,snapshot_number,group_text_prefix);
      strcpy(fp_out->filename_properties_base,fp_out->filename_properties_root);
      strcpy(fp_out->filename_profiles_base,  fp_out->filename_profiles_root);
      strip_path(fp_out->filename_properties_base);
      strip_path(fp_out->filename_profiles_base);

      // Try reading a multifile first ...
      sprintf(filename_properties,"%s/%s.%d",fp_out->filename_properties_root,fp_out->filename_properties_base,0);
      fp_out->fp_properties=fopen(filename_properties,"r");
      if(fp_out->fp_properties==NULL){
         sprintf(filename_properties,"%s",fp_out->filename_properties_root);
         // ... if we didn't find a multi-file, try reading a single file ...
         fp_out->fp_properties=fopen(filename_properties,"r");
         if(fp_out->fp_properties==NULL){
           r_val=TRUE;
           SID_trap_error("Could not open catalog {%s} snapshot #%03d {%s}.",ERROR_IO_OPEN,filename_catalog_root,snapshot_number,filename_properties);
         }
         // ... we found a single file.  Set flags.
         else
           fp_out->flag_multifile=FALSE;
      }
      // ... we found a multi-file.  Set flags.
      else
         fp_out->flag_multifile=TRUE;

      // Load/set header information
      if(fp_out->fp_properties!=NULL){
         fread(&(fp_out->i_file),       sizeof(int),1,fp_out->fp_properties);
         fread(&(fp_out->n_files),      sizeof(int),1,fp_out->fp_properties);
         fread(&(fp_out->n_halos_file), sizeof(int),1,fp_out->fp_properties);
         fread(&(fp_out->n_halos_total),sizeof(int),1,fp_out->fp_properties);
         fclose(fp_out->fp_properties);
         fp_out->fp_properties=NULL;
      }
   }
   SID_Bcast(fp_out,sizeof(fp_catalog_info),MASTER_RANK,SID.COMM_WORLD);

   // Initialize things by opening the first file.  Even if we want a halo
   //   that's deep in the list, we have to scan all the headers (starting
   //   with the first) to find where it is.
   int r_val2;
   fp_out->fp_properties=NULL;
   fp_out->fp_profiles  =NULL;
   fp_out->i_file       =0;
   fp_out->i_halo       =0;
   if((r_val2=fopen_nth_catalog_file(fp_out,0)))
      SID_trap_error("Error opening catalog file.",ERROR_IO_OPEN);
   if(r_val2>r_val)
      r_val=r_val2;

   return(r_val);
}

