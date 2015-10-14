#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpSims_Swap_Endian.h>

int swap_endian_smooth(const char *filename_in_root,const char *filename_out_root,int region_number,int snap_number,int mode,int IDs_byte_size){

  if(region_number<0)
     SID_log("Swapping endian of full snapshot smooth file...",SID_LOG_OPEN);
  else
     SID_log("Swapping endian of region #%03d smooth file...",SID_LOG_OPEN,region_number);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_catalogs_properties_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  if(region_number<0){
     sprintf(filename_in, "%s/smooth/smooth_%03d/smooth_%03d.%d",filename_in_root, snap_number,snap_number,0);
     sprintf(filename_out,"%s/smooth/smooth_%03d/smooth_%03d.%d",filename_out_root,snap_number,snap_number,0);
  }
  else{
     sprintf(filename_in, "%s/smooth/smooth_region%03d_%03d/smooth_region%03d_%03d.%d",filename_in_root, region_number,snap_number,region_number,snap_number,0);
     sprintf(filename_out,"%s/smooth/smooth_region%03d_%03d/smooth_region%03d_%03d.%d",filename_out_root,region_number,snap_number,region_number,snap_number,0);
  }

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  int   flag_type=FALSE;
  if((fp_in=fopen(filename_in,"r"))==NULL){
     if(region_number<0)
        sprintf(filename_in, "%s/smooth/smooth_%03d",filename_in_root,snap_number);
     else
        sprintf(filename_in, "%s/smooth/smooth_region%03d_%03d",filename_in_root,region_number,snap_number);
     if((fp_in=fopen(filename_in,"r"))==NULL){
        SID_log("not present.",SID_LOG_CLOSE);
        return(FALSE);
     }
     flag_type=TRUE;
  }
  else{
     char filename_dir[MAX_FILENAME_LENGTH];
     if(region_number<0)
        sprintf(filename_dir,"%s/smooth/smooth_%03d",filename_out_root,snap_number);
     else
        sprintf(filename_dir,"%s/smooth/smooth_region%03d_%03d",filename_out_root,region_number,snap_number);
     mkdir(filename_dir,02755);
  }

  // Read the needed header information
  int       n_particles_file;
  int       offset;
  long long n_particles_total;
  int       n_files;
  fread(&n_particles_file, sizeof(int),      1,fp_in);
  fread(&offset,           sizeof(int),      1,fp_in);
  fread(&n_particles_total,sizeof(long long),1,fp_in);
  fread(&n_files,          sizeof(int),      1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
     swap_endian((char *)(&n_particles_file), 1,sizeof(int));
     swap_endian((char *)(&offset),           1,sizeof(int));
     swap_endian((char *)(&n_particles_total),1,sizeof(long long));
     swap_endian((char *)(&n_files),          1,sizeof(int));
  }
  fclose(fp_in);

  // Sanity check
  if(flag_type && n_files!=1) SID_trap_error("Single file format dataset {%s} has n_files=%d (ie !=1) in its header.",ERROR_LOGIC,filename_in,n_files);

  int i_file;
  for(i_file=0;i_file<n_files;i_file++){
     // Open the files
     if(flag_type){
        if(region_number<0){
           sprintf(filename_in, "%s/smooth/smooth_%03d",filename_in_root, snap_number);
           sprintf(filename_out,"%s/smooth/smooth_%03d",filename_out_root,snap_number);
        }
        else{
           sprintf(filename_in, "%s/smooth/smooth_region%03d_%03d",filename_in_root, region_number,snap_number);
           sprintf(filename_out,"%s/smooth/smooth_region%03d_%03d",filename_out_root,region_number,snap_number);
        }
     }
     else{
        if(region_number<0){
           sprintf(filename_in, "%s/smooth/smooth_%03d/smooth_%03d.%d",filename_in_root, snap_number,snap_number,i_file);
           sprintf(filename_out,"%s/smooth/smooth_%03d/smooth_%03d.%d",filename_out_root,snap_number,snap_number,i_file);
        }
        else{
           sprintf(filename_in, "%s/smooth/smooth_region%03d_%03d/smooth_region%03d_%03d.%d",filename_in_root, 
                   region_number,snap_number,region_number,snap_number,i_file);
           sprintf(filename_out,"%s/smooth/smooth_region%03d_%03d/smooth_region%03d_%03d.%d",filename_out_root,
                   region_number,snap_number,region_number,snap_number,i_file);
        }
     }
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
     if((fp_out=fopen(filename_out,"w"))==NULL)
        SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

     // Process the header
     fread(&n_particles_file, sizeof(int),      1,fp_in);
     fread(&offset,           sizeof(int),      1,fp_in);
     fread(&n_particles_total,sizeof(long long),1,fp_in);
     fread(&n_files,          sizeof(int),      1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
        swap_endian((char *)(&n_particles_file), 1,sizeof(int));
        swap_endian((char *)(&offset),           1,sizeof(int));
        swap_endian((char *)(&n_particles_total),1,sizeof(long long));
        swap_endian((char *)(&n_files),          1,sizeof(int));
     }
     rewind(fp_in);
     rewrite_swap_endian(fp_in,fp_out,2,sizeof(int),      NULL);
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(long long),NULL);
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),      NULL);

     // Create buffer
     char *buffer=(char *)SID_malloc(32);

     // Process smoothing lengths
     for(int i_particle=0;i_particle<n_particles_file;i_particle++)
        rewrite_swap_endian(fp_in,fp_out,1,sizeof(float),buffer);

     // Process densities
     for(int i_particle=0;i_particle<n_particles_file;i_particle++)
        rewrite_swap_endian(fp_in,fp_out,1,sizeof(float),buffer);

     // Process velocity dispersions
     for(int i_particle=0;i_particle<n_particles_file;i_particle++)
        rewrite_swap_endian(fp_in,fp_out,1,sizeof(float),buffer);

     // Process IDs
     for(int i_particle=0;i_particle<n_particles_file;i_particle++)
        rewrite_swap_endian(fp_in,fp_out,1,IDs_byte_size,buffer);

     // Free buffer
     SID_free(SID_FARG buffer);

     // Close files
     fclose(fp_in);
     fclose(fp_out);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);
  return(TRUE);
}

