#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpSwapEndian.h>

void swap_endian_catalogs_properties_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,const char *prefix,int snap_number,int mode);
void swap_endian_catalogs_properties_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,const char *prefix,int snap_number,int mode){

  SID_log("Swapping endian of %sgroup properties...",SID_LOG_OPEN,prefix);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_catalogs_properties_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_properties/%s_%03d.catalog_%sgroups_properties.%d",filename_in_root, filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);
  sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_properties/%s_%03d.catalog_%sgroups_properties.%d",filename_out_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_properties",filename_in_root,filename_halo_type,snap_number,prefix);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in_root);
     flag_type=TRUE;
  }
  else{
     char filename_dir[MAX_FILENAME_LENGTH];
     sprintf(filename_dir,"%s/%s_%03d.catalog_%sgroups_properties",filename_out_root,filename_halo_type,snap_number,prefix);
     mkdir(filename_dir,02755);
  }

  // Read the needed header information
  int i_file_in;
  int n_files;
  fread_verify(&i_file_in,sizeof(int),1,fp_in);
  fread_verify(&n_files,  sizeof(int),1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
     swap_endian((char *)(&n_files),1,sizeof(int));
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  int i_file;
  for(i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type){
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_properties",filename_in_root, filename_halo_type,snap_number,prefix);
        sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_properties",filename_out_root,filename_halo_type,snap_number,prefix);
     }
     else{
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_properties/%s_%03d.catalog_%sgroups_properties.%d",filename_in_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
        sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_properties/%s_%03d.catalog_%sgroups_properties.%d",filename_out_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
     }
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
     if((fp_out=fopen(filename_out,"w"))==NULL)
        SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

     // Read the header and rewind
     int n_halos_file;
     int n_halos_total;
     fread_verify(&i_file_in,    sizeof(int),1,fp_in);
     fread_verify(&n_files,      sizeof(int),1,fp_in);
     fread_verify(&n_halos_file, sizeof(int),1,fp_in);
     fread_verify(&n_halos_total,sizeof(int),1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
        swap_endian((char *)(&i_file_in),    1,sizeof(int));
        swap_endian((char *)(&n_files),      1,sizeof(int));
        swap_endian((char *)(&n_halos_file), 1,sizeof(int));
        swap_endian((char *)(&n_halos_total),1,sizeof(int));
     }
     rewind(fp_in);

     // Process the file
     rewrite_swap_endian(fp_in,fp_out,4,sizeof(int),NULL);
     for(int i_halo=0;i_halo<n_halos_file;i_halo++){
        halo_properties_info properties;
        fread_verify(&properties,sizeof(halo_properties_info),1,fp_in);
        swap_endian((char *)(&(properties.id_MBP)),             1,sizeof(long long));
        swap_endian((char *)(&(properties.M_vir)),              1,sizeof(double));
        swap_endian((char *)(&(properties.n_particles)),        1,sizeof(int));
        swap_endian((char *)(&(properties.position_COM)),       3,sizeof(float));
        swap_endian((char *)(&(properties.position_MBP)),       3,sizeof(float));
        swap_endian((char *)(&(properties.velocity_COM)),       3,sizeof(float));
        swap_endian((char *)(&(properties.velocity_MBP)),       3,sizeof(float));
        swap_endian((char *)(&(properties.R_vir)),              1,sizeof(float));
        swap_endian((char *)(&(properties.R_halo)),             1,sizeof(float));
        swap_endian((char *)(&(properties.R_max)),              1,sizeof(float));
        swap_endian((char *)(&(properties.V_max)),              1,sizeof(float));
        swap_endian((char *)(&(properties.sigma_v)),            1,sizeof(float));
        swap_endian((char *)(&(properties.spin)),               3,sizeof(float));
        swap_endian((char *)(&(properties.q_triaxial)),         1,sizeof(float));
        swap_endian((char *)(&(properties.s_triaxial)),         1,sizeof(float));
        swap_endian((char *)(&(properties.shape_eigen_vectors)),9,sizeof(float));
        fwrite(&properties,sizeof(halo_properties_info),1,fp_out);
     }

     // Close files
     fclose(fp_in);
     fclose(fp_out);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);
}

void swap_endian_catalogs_profiles_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,const char *prefix,int snap_number,int mode);
void swap_endian_catalogs_profiles_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,const char *prefix,int snap_number,int mode){

  SID_log("Swapping endian of %sgroup profiles...",SID_LOG_OPEN,prefix);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_catalogs_profiles_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_profiles/%s_%03d.catalog_%sgroups_profiles.%d",filename_in_root, filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);
  sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_profiles/%s_%03d.catalog_%sgroups_profiles.%d",filename_out_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_profiles",filename_in_root,filename_halo_type,snap_number,prefix);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in_root);
     flag_type=TRUE;
  }
  else{
     char filename_dir[MAX_FILENAME_LENGTH];
     sprintf(filename_dir,"%s/%s_%03d.catalog_%sgroups_profiles",filename_out_root,filename_halo_type,snap_number,prefix);
     mkdir(filename_dir,02755);
  }

  // Read the needed header information
  int i_file_in;
  int n_files;
  fread_verify(&i_file_in,sizeof(int),1,fp_in);
  fread_verify(&n_files,  sizeof(int),1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
     swap_endian((char *)(&n_files),1,sizeof(int));
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  for(int i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type){
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_profiles",filename_in_root, filename_halo_type,snap_number,prefix);
        sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_profiles",filename_out_root,filename_halo_type,snap_number,prefix);
     }
     else{
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_profiles/%s_%03d.catalog_%sgroups_profiles.%d",filename_in_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
        sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_profiles/%s_%03d.catalog_%sgroups_profiles.%d",filename_out_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
     }
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
     if((fp_out=fopen(filename_out,"w"))==NULL)
        SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

     // Read the header and rewind
     int n_halos_file;
     int n_halos_total;
     fread_verify(&i_file_in,    sizeof(int),1,fp_in);
     fread_verify(&n_files,      sizeof(int),1,fp_in);
     fread_verify(&n_halos_file, sizeof(int),1,fp_in);
     fread_verify(&n_halos_total,sizeof(int),1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
        swap_endian((char *)(&i_file_in),    1,sizeof(int));
        swap_endian((char *)(&n_files),      1,sizeof(int));
        swap_endian((char *)(&n_halos_file), 1,sizeof(int));
        swap_endian((char *)(&n_halos_total),1,sizeof(int));
     }
     rewind(fp_in);

     // Process the file
     rewrite_swap_endian(fp_in,fp_out,4,sizeof(int),NULL);
     for(int i_halo=0;i_halo<n_halos_file;i_halo++){
        int n_bins;
        int n_bins_in;
        int n_bins_out;
        fread_verify(&n_bins_in,sizeof(int),1,fp_in);
        n_bins    =n_bins_in;
        n_bins_out=n_bins_in;
        swap_endian((char *)(&n_bins_out),1,sizeof(int));
        if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
           swap_endian((char *)(&n_bins),1,sizeof(int));
        fwrite(&n_bins_out,sizeof(int),1,fp_out);
        for(int i_bin=0;i_bin<n_bins;i_bin++){
           halo_profile_bin_info profile_bin;
           fread_verify(&profile_bin,sizeof(halo_profile_bin_info),1,fp_in);
           swap_endian((char *)(&(profile_bin.r_med)),              1,sizeof(float));
           swap_endian((char *)(&(profile_bin.r_max)),              1,sizeof(float));
           swap_endian((char *)(&(profile_bin.n_particles)),        1,sizeof(int));
           swap_endian((char *)(&(profile_bin.rho)),                1,sizeof(float));
           swap_endian((char *)(&(profile_bin.M_r)),                1,sizeof(double));
           swap_endian((char *)(&(profile_bin.overdensity)),        1,sizeof(float));
           swap_endian((char *)(&(profile_bin.position_COM)),       3,sizeof(float));
           swap_endian((char *)(&(profile_bin.velocity_COM)),       3,sizeof(float));
           swap_endian((char *)(&(profile_bin.sigma_rad)),          1,sizeof(float));
           swap_endian((char *)(&(profile_bin.sigma_tan)),          1,sizeof(float));
           swap_endian((char *)(&(profile_bin.sigma_tot)),          1,sizeof(float));
           swap_endian((char *)(&(profile_bin.spin)),               3,sizeof(float));
           swap_endian((char *)(&(profile_bin.q_triaxial)),         1,sizeof(float));
           swap_endian((char *)(&(profile_bin.s_triaxial)),         1,sizeof(float));
           swap_endian((char *)(&(profile_bin.shape_eigen_vectors)),9,sizeof(float));
           fwrite(&profile_bin,sizeof(halo_profile_bin_info),1,fp_out);
        }
     }

     // Close files
     fclose(fp_in);
     fclose(fp_out);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);
}

int swap_endian_catalogs_SO_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,const char *prefix,int snap_number,int mode);
int swap_endian_catalogs_SO_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,const char *prefix,int snap_number,int mode){
  SID_log("Swapping endian of %sgroup SO properties...",SID_LOG_OPEN,prefix);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_catalogs_properties_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_SO/%s_%03d.catalog_%sgroups_SO.%d",filename_in_root, filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);
  sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_SO/%s_%03d.catalog_%sgroups_SO.%d",filename_out_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,0);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     sprintf(filename_in,"%s/%s_%03d.catalog_%sgroups_SO",filename_in_root,filename_halo_type,snap_number,prefix);
     if((fp_in=fopen(filename_in,"r"))==NULL){
        SID_log("not present.",SID_LOG_CLOSE);
        return(FALSE);
     }
     flag_type=TRUE;
  }
  else{
     char filename_dir[MAX_FILENAME_LENGTH];
     sprintf(filename_dir,"%s/%s_%03d.catalog_%sgroups_SO",filename_out_root,filename_halo_type,snap_number,prefix);
     mkdir(filename_dir,02755);
  }

  // Read the needed header information
  int i_file_in;
  int n_files;
  fread_verify(&i_file_in,sizeof(int),1,fp_in);
  fread_verify(&n_files,  sizeof(int),1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
     swap_endian((char *)(&n_files),1,sizeof(int));
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  // Process each file in turn
  for(int i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type){
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_SO",filename_in_root, filename_halo_type,snap_number,prefix);
        sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_SO",filename_out_root,filename_halo_type,snap_number,prefix);
     }
     else{
        sprintf(filename_in, "%s/%s_%03d.catalog_%sgroups_SO/%s_%03d.catalog_%sgroups_SO.%d",filename_in_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
        sprintf(filename_out,"%s/%s_%03d.catalog_%sgroups_SO/%s_%03d.catalog_%sgroups_SO.%d",filename_out_root,filename_halo_type,snap_number,prefix,filename_halo_type,snap_number,prefix,i_file);
     }
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
     if((fp_out=fopen(filename_out,"w"))==NULL)
        SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

     // Read the header and rewind
     int n_halos_file;
     int n_halos_total;
     fread_verify(&i_file_in,    sizeof(int),1,fp_in);
     fread_verify(&n_files,      sizeof(int),1,fp_in);
     fread_verify(&n_halos_file, sizeof(int),1,fp_in);
     fread_verify(&n_halos_total,sizeof(int),1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
        swap_endian((char *)(&i_file_in),    1,sizeof(int));
        swap_endian((char *)(&n_files),      1,sizeof(int));
        swap_endian((char *)(&n_halos_file), 1,sizeof(int));
        swap_endian((char *)(&n_halos_total),1,sizeof(int));
     }
     rewind(fp_in);

     // Process the file
     rewrite_swap_endian(fp_in,fp_out,4,sizeof(int),NULL);
     for(int i_halo=0;i_halo<(6*n_halos_file);i_halo++)
        rewrite_swap_endian(fp_in,fp_out,1,sizeof(float),NULL);

     // Close files
     fclose(fp_in);
     fclose(fp_out);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);

}

void swap_endian_catalogs(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode){
   for(int i_type=0;i_type<2;i_type++){
      char prefix[8];
      if(i_type==0)
         sprintf(prefix,"");
      else
         sprintf(prefix,"sub");
      swap_endian_catalogs_properties_local(filename_in_root,filename_out_root,filename_halo_type,prefix,snap_number,mode);
      swap_endian_catalogs_profiles_local  (filename_in_root,filename_out_root,filename_halo_type,prefix,snap_number,mode);
      if(i_type==0)
         swap_endian_catalogs_SO_local(filename_in_root,filename_out_root,filename_halo_type,prefix,snap_number,mode);
   }
}

