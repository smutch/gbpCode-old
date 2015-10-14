#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpSims_Swap_Endian.h>

void check_integrity_catalogs_properties_local(const char *filename_in_root,const char *filename_halo_type,const char *prefix,int snap_number);
void check_integrity_catalogs_properties_local(const char *filename_in_root,const char *filename_halo_type,const char *prefix,int snap_number){

  SID_log("Swapping endian of %sgroup properties...",SID_LOG_OPEN,prefix);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_properties/%s_%03d.catalog_%sgroups_properties.%d",filename_in_root, filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);

  // Open input files
  FILE *fp_in =NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_properties",filename_in_root,filename_halo_type,snap_number,prefix);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in_root);
     flag_type=TRUE;
  }

  // Read the needed header information
  int i_file_in;
  int n_files;
  fread_verify(&i_file_in,sizeof(int),1,fp_in);
  fread_verify(&n_files,  sizeof(int),1,fp_in);
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  int i_file;
  for(i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type)
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_properties",filename_in_root, filename_halo_type,snap_number,prefix);
     else
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_properties/%s_%03d.catalog_%sgroups_properties.%d",filename_in_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

     // Read the header and rewind
     int n_halos_file;
     int n_halos_total;
     fread_verify(&i_file_in,    sizeof(int),1,fp_in);
     fread_verify(&n_files,      sizeof(int),1,fp_in);
     fread_verify(&n_halos_file, sizeof(int),1,fp_in);
     fread_verify(&n_halos_total,sizeof(int),1,fp_in);

     // Process the file
     for(int i_halo=0;i_halo<n_halos_file;i_halo++){
        halo_properties_info properties;
        fread_verify(&properties,sizeof(halo_properties_info),1,fp_in);
     }

     // Check that we are at the end of the file
     char test;fread(&test,1,1,fp_in);
     if(!feof(fp_in))
        SID_trap_error("There are stray bytes at the end of {%s}.",ERROR_LOGIC,filename_in);

     // Close files
     fclose(fp_in);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);
}

void check_integrity_catalogs_profiles_local(const char *filename_in_root,const char *filename_halo_type,const char *prefix,int snap_number);
void check_integrity_catalogs_profiles_local(const char *filename_in_root,const char *filename_halo_type,const char *prefix,int snap_number){

  SID_log("Swapping endian of %sgroup profiles...",SID_LOG_OPEN,prefix);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_profiles/%s_%03d.catalog_%sgroups_profiles.%d",filename_in_root, filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);

  // Open input files
  FILE *fp_in =NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_profiles",filename_in_root,filename_halo_type,snap_number,prefix);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in_root);
     flag_type=TRUE;
  }

  // Read the needed header information
  int i_file_in;
  int n_files;
  fread_verify(&i_file_in,sizeof(int),1,fp_in);
  fread_verify(&n_files,  sizeof(int),1,fp_in);
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  for(int i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type)
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_profiles",filename_in_root, filename_halo_type,snap_number,prefix);
     else
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_profiles/%s_%03d.catalog_%sgroups_profiles.%d",filename_in_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

     // Read the header and rewind
     int n_halos_file;
     int n_halos_total;
     fread_verify(&i_file_in,    sizeof(int),1,fp_in);
     fread_verify(&n_files,      sizeof(int),1,fp_in);
     fread_verify(&n_halos_file, sizeof(int),1,fp_in);
     fread_verify(&n_halos_total,sizeof(int),1,fp_in);

     // Process the file
     for(int i_halo=0;i_halo<n_halos_file;i_halo++){
        int n_bins;
        fread_verify(&n_bins,sizeof(int),1,fp_in);
        for(int i_bin=0;i_bin<n_bins;i_bin++){
           halo_profile_bin_info profile_bin;
           fread_verify(&profile_bin,sizeof(halo_profile_bin_info),1,fp_in);
        }
     }

     // Check that we are at the end of the file
     char test;fread(&test,1,1,fp_in);
     if(!feof(fp_in))
        SID_trap_error("There are stray bytes at the end of {%s}.",ERROR_LOGIC,filename_in);

     // Close files
     fclose(fp_in);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);
}

int check_integrity_catalogs_SO_local(const char *filename_in_root,const char *filename_halo_type,const char *prefix,int snap_number);
int check_integrity_catalogs_SO_local(const char *filename_in_root,const char *filename_halo_type,const char *prefix,int snap_number){
  SID_log("Swapping endian of %sgroup SO properties...",SID_LOG_OPEN,prefix);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_SO/%s_%03d.catalog_%sgroups_SO.%d",filename_in_root, filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);

  // Open input files
  FILE *fp_in =NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_SO",filename_in_root,filename_halo_type,snap_number,prefix);
     if((fp_in=fopen(filename_in,"r"))==NULL){
        SID_log("not present.",SID_LOG_CLOSE);
        return(FALSE);
     }
     flag_type=TRUE;
  }

  // Read the needed header information
  int i_file_in;
  int n_files;
  fread_verify(&i_file_in,sizeof(int),1,fp_in);
  fread_verify(&n_files,  sizeof(int),1,fp_in);
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  // Process each file in turn
  for(int i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type)
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_SO",filename_in_root, filename_halo_type,snap_number,prefix);
     else
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_SO/%s_%03d.catalog_%sgroups_SO.%d",filename_in_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

     // Read the header and rewind
     int n_halos_file;
     int n_halos_total;
     fread_verify(&i_file_in,    sizeof(int),1,fp_in);
     fread_verify(&n_files,      sizeof(int),1,fp_in);
     fread_verify(&n_halos_file, sizeof(int),1,fp_in);
     fread_verify(&n_halos_total,sizeof(int),1,fp_in);

     // Process the file
     for(int i_halo=0;i_halo<(6*n_halos_file);i_halo++){
        float buffer;
        fread_verify(&buffer,sizeof(float),1,fp_in);
     }

     // Check that we are at the end of the file
     char test;fread(&test,1,1,fp_in);
     if(!feof(fp_in))
        SID_trap_error("There are stray bytes at the end of {%s}.",ERROR_LOGIC,filename_in);

     // Close files
     fclose(fp_in);
  }
  SID_log("Done.",SID_LOG_CLOSE);

}

void check_integrity_catalogs(const char *filename_in_root,const char *filename_halo_type,int snap_number){
   for(int i_type=0;i_type<2;i_type++){
      char prefix[8];
      if(i_type==0)
         sprintf(prefix,"");
      else
         sprintf(prefix,"sub");
      check_integrity_catalogs_properties_local(filename_in_root,filename_halo_type,prefix,snap_number);
      check_integrity_catalogs_profiles_local  (filename_in_root,filename_halo_type,prefix,snap_number);
      if(i_type==0)
         check_integrity_catalogs_SO_local(filename_in_root,filename_halo_type,prefix,snap_number);
   }
}

