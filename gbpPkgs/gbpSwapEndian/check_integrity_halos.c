#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSwapEndian.h>

void check_integrity_halos_groups_local(const char *filename_in_root,const char *filename_halo_type,int snap_number);
void check_integrity_halos_groups_local(const char *filename_in_root,const char *filename_halo_type,int snap_number){

  SID_log("Checking group file...",SID_LOG_OPEN);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_groups",filename_in_root, filename_halo_type,snap_number);

  // Open input files
  FILE *fp_in =NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

  // Read the needed header information and rewind
  int n_groups;
  int offset_size_bytes;
  fread_verify(&n_groups,         sizeof(int),1,fp_in);
  fread_verify(&offset_size_bytes,sizeof(int),1,fp_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*offset_size_bytes);

  // Process the file
  for(int i_group=0;i_group<n_groups;i_group++)
     fread_verify(buffer,sizeof(int),1,fp_in);
  for(int i_group=0;i_group<n_groups;i_group++)
     fread_verify(buffer,offset_size_bytes,1,fp_in);
  for(int i_group=0;i_group<n_groups;i_group++)
     fread_verify(buffer,sizeof(int),1,fp_in);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Check that we are at the end of the file
  char test;fread(&test,1,1,fp_in);
  if(!feof(fp_in))
     SID_trap_error("There are stray bytes at the end of {%s}.",ERROR_LOGIC,filename_in);

  // Close files
  fclose(fp_in);
  
  SID_log("Done.",SID_LOG_CLOSE);
}

void check_integrity_halos_subgroups_local(const char *filename_in_root,const char *filename_halo_type,int snap_number);
void check_integrity_halos_subgroups_local(const char *filename_in_root,const char *filename_halo_type,int snap_number){

  SID_log("Swapping endian of subgroup file...",SID_LOG_OPEN);

  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_subgroups",filename_in_root, filename_halo_type,snap_number);

  // Open input files
  FILE *fp_in =NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

  // Read the needed header information and rewind
  int n_subgroups;
  int offset_size_bytes;
  fread_verify(&n_subgroups,      sizeof(int),1,fp_in);
  fread_verify(&offset_size_bytes,sizeof(int),1,fp_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*offset_size_bytes);

  // Process the file
  for(int i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++)
     fread_verify(buffer,sizeof(int),1,fp_in);
  for(int i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++)
     fread_verify(buffer,offset_size_bytes,1,fp_in);
  for(int i_subgroup=0;i_subgroup<n_subgroups;i_subgroup++)
     fread_verify(buffer,sizeof(int),1,fp_in);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Check that we are at the end of the file
  char test;fread(&test,1,1,fp_in);
  if(!feof(fp_in))
     SID_trap_error("There are stray bytes at the end of {%s}.",ERROR_LOGIC,filename_in);

  // Close files
  fclose(fp_in);

  SID_log("Done.",SID_LOG_CLOSE);
}

void check_integrity_halos_particles_local(const char *filename_in_root,const char *filename_halo_type,int snap_number);
void check_integrity_halos_particles_local(const char *filename_in_root,const char *filename_halo_type,int snap_number){

  SID_log("Swapping endian of particles file...",SID_LOG_OPEN);
   
  // Set filenames
  char filename_in[MAX_FILENAME_LENGTH];
  sprintf(filename_in, "%s/%s_%03d.catalog_particles",filename_in_root, filename_halo_type,snap_number);

  // Open input files
  FILE *fp_in =NULL;
  if((fp_in=fopen(filename_in,"r"))==NULL)
     SID_trap_error("Could not open {%s} for reading.",ERROR_IO_OPEN,filename_in);

  // Read the needed header information and rewind
  int       particle_size_bytes;
  int       n_particles_int;
  long long n_particles_long_long;
  size_t    n_particles;
  fread_verify(&particle_size_bytes,sizeof(int),1,fp_in);
  if(particle_size_bytes==sizeof(int)){
     fread_verify(&n_particles_int,sizeof(int),1,fp_in);
     n_particles=(size_t)n_particles_int;
  }
  else if(particle_size_bytes==sizeof(long long)){
     fread_verify(&n_particles_long_long,sizeof(long long),1,fp_in);
     n_particles=(size_t)n_particles_long_long;
  }
  else
     SID_trap_error("Invalid particle ID size (%d) in {%s}.",ERROR_LOGIC,particle_size_bytes,filename_in);

  // Create a read buffer
  char *buffer=(char *)SID_malloc(sizeof(char)*particle_size_bytes);

  // Process the file
  for(size_t i_particle=0;i_particle<n_particles;i_particle++)
     fread_verify(buffer,particle_size_bytes,1,fp_in);

  // Free the read buffer
  SID_free(SID_FARG buffer);

  // Check that we are at the end of the file
  char test;fread(&test,1,1,fp_in);
  if(!feof(fp_in))
     SID_trap_error("There are stray bytes at the end of {%s}.",ERROR_LOGIC,filename_in);

  // Close files
  fclose(fp_in);

  SID_log("Done.",SID_LOG_CLOSE);
}

void check_integrity_halos(const char *filename_in_root,const char *filename_halo_type,int snap_number){
   check_integrity_halos_groups_local   (filename_in_root,filename_halo_type,snap_number);
   check_integrity_halos_subgroups_local(filename_in_root,filename_halo_type,snap_number);
   check_integrity_halos_particles_local(filename_in_root,filename_halo_type,snap_number);
}

