#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSims_Swap_Endian.h>

void swap_endian_halos_groups_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode);
void swap_endian_halos_groups_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode){

  SID_log("Swapping endian of group file...",SID_LOG_OPEN);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_halos_groups_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_groups",filename_in_root, filename_halo_type,snap_number);
  sprintf(filename_out,"%s/%s_%03d.catalog_groups",filename_out_root,filename_halo_type,snap_number);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
  if((fp_out=fopen(filename_out,"w"))==NULL)
     SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

  // Read the needed header information and rewind
  int n_groups;
  int offset_size_bytes;
  fread(&n_groups,         sizeof(int),1,fp_in);
  fread(&offset_size_bytes,sizeof(int),1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
     swap_endian((char *)(&n_groups),         1,sizeof(int));
     swap_endian((char *)(&offset_size_bytes),1,sizeof(int));
  }
  rewind(fp_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*offset_size_bytes);

  // Process the file
  rewrite_swap_endian(fp_in,fp_out,2,sizeof(int),   buffer);
  for(int i_group=0;i_group<n_groups;i_group++)
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);
  for(int i_group=0;i_group<n_groups;i_group++)
     rewrite_swap_endian(fp_in,fp_out,1,offset_size_bytes,buffer);
  for(int i_group=0;i_group<n_groups;i_group++)
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Close files
  fclose(fp_in);
  fclose(fp_out);
  
  SID_log("Done.",SID_LOG_CLOSE);
}

void swap_endian_halos_subgroups_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode);
void swap_endian_halos_subgroups_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode){

  SID_log("Swapping endian of subgroup file...",SID_LOG_OPEN);

  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_halos_subgroups_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_subgroups",filename_in_root, filename_halo_type,snap_number);
  sprintf(filename_out,"%s/%s_%03d.catalog_subgroups",filename_out_root,filename_halo_type,snap_number);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
  if((fp_out=fopen(filename_out,"w"))==NULL)
     SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

  // Read the needed header information and rewind
  int n_subgroups;
  int offset_size_bytes;
  fread(&n_subgroups,      sizeof(int),1,fp_in);
  fread(&offset_size_bytes,sizeof(int),1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE)){
     swap_endian((char *)(&n_subgroups),      1,sizeof(int));
     swap_endian((char *)(&offset_size_bytes),1,sizeof(int));
  }
  rewind(fp_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*offset_size_bytes);

  // Process the file
  rewrite_swap_endian(fp_in,fp_out,2,sizeof(int),   buffer);
  for(int i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++)
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);
  for(int i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++)
     rewrite_swap_endian(fp_in,fp_out,1,offset_size_bytes,buffer);
  for(int i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++)
     rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),buffer);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Close files
  fclose(fp_in);
  fclose(fp_out);

  SID_log("Done.",SID_LOG_CLOSE);
}

void swap_endian_halos_particles_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode);
void swap_endian_halos_particles_local(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode){

  SID_log("Swapping endian of particles file...",SID_LOG_OPEN);
   
  // Sanity check
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE) && check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_FROM_NATIVE))
     SID_trap_error("Invalid mode flag (%d) in swap_endian_halos_particles_local().",ERROR_LOGIC,mode);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_particles",filename_in_root, filename_halo_type,snap_number);
  sprintf(filename_out,"%s/%s_%03d.catalog_particles",filename_out_root,filename_halo_type,snap_number);

  // Open input and output files
  FILE *fp_in =NULL;
  FILE *fp_out=NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);
  if((fp_out=fopen(filename_out,"w"))==NULL)
     SID_trap_error("Could not open {%s} for writing.",ERROR_IO_OPEN,filename_out);

  // Read the needed header information and rewind
  int       particle_size_bytes;
  int       n_particles_int;
  long long n_particles_long_long;
  size_t    n_particles;
  fread(&particle_size_bytes,sizeof(int),1,fp_in);
  if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
     swap_endian((char *)(&particle_size_bytes),1,sizeof(int));
  if(particle_size_bytes==sizeof(int)){
     fread(&n_particles_int,sizeof(int),1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
        swap_endian((char *)(&n_particles_int),1,sizeof(int));
     n_particles=(size_t)n_particles_int;
  }
  else if(particle_size_bytes==sizeof(long long)){
     fread(&n_particles_long_long,sizeof(long long),1,fp_in);
     if(check_mode_for_flag(mode,SWAP_SSIMPL_ENDIAN_TO_NATIVE))
        swap_endian((char *)(&n_particles_long_long),1,sizeof(long long));
     n_particles=(size_t)n_particles_long_long;
  }
  else
     SID_trap_error("Invalid particle ID size (%d) in {%s}.",ERROR_LOGIC,particle_size_bytes,filename_in);
  rewind(fp_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*particle_size_bytes);

  // Process the file
  rewrite_swap_endian(fp_in,fp_out,1,sizeof(int),         buffer);
  rewrite_swap_endian(fp_in,fp_out,1,particle_size_bytes,buffer);
  for(size_t i_particle=0;i_particle<n_particles;i_particle++)
     rewrite_swap_endian(fp_in,fp_out,1,particle_size_bytes,buffer);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Close files
  fclose(fp_in);
  fclose(fp_out);

  SID_log("Done.",SID_LOG_CLOSE);
}

void swap_endian_halos(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode){
   swap_endian_halos_groups_local   (filename_in_root,filename_out_root,filename_halo_type,snap_number,mode);
   swap_endian_halos_subgroups_local(filename_in_root,filename_out_root,filename_halo_type,snap_number,mode);
   swap_endian_halos_particles_local(filename_in_root,filename_out_root,filename_halo_type,snap_number,mode);
}

