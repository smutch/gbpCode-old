#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpHalos.h>

void read_groups_AHF(char        *filename_groups_root,
                     int          i_file,
                     int          mode,
                     plist_info  *plist,
                     char        *catalog_name){
  char     filename_ids[256];
  char     filename_number[256];
  FILE    *fp;
  FILE    *fp_2;
  int      n_groups;
  size_t   i_particle;
  size_t   j_particle;
  size_t   k_particle;
  int     *n_particles_group;
  int     *group_offsets;
  int      n_particles_inemp;
  big_int  n_particles_in;
  big_int *particle_ids;
  big_int  n_particles;
  big_int  dummy;
  big_int  n_groups_file;
  size_t   i_group;
  size_t   j_group;

  // Set AHF file name
  if(i_file<10)
    sprintf(filename_number,"00%1d",i_file);
  else if(i_file<100)
    sprintf(filename_number,"0%2d", i_file);
  else
    sprintf(filename_number,"%3d", i_file);
  sprintf(filename_ids,"%s_%s.AHF_particles",filename_groups_root,filename_number);

  // Open group particle list file
  if((fp=fopen(filename_ids,"r"))!=NULL){

    // Read group particle list file
    SID_log("Reading particle file {%s}...",SID_LOG_OPEN,filename_ids);

    // Count particles and number of groups
    n_particles=0;
    n_groups   =0;
    fscanf(fp,"%lld",&n_groups_file);
    while(!feof(fp)){
      for(j_group=0;j_group<n_groups_file;j_group++){
        fscanf(fp,"%lld",&n_particles_in);
        n_groups++;
        n_particles+=n_particles_in;
        for(j_particle=0;j_particle<n_particles_in;j_particle++){
          fscanf(fp,"%lld",&dummy);
        }
      }
      fscanf(fp,"%lld",&n_groups_file);
    }
    SID_log("n_groups   =%d",  SID_LOG_COMMENT,n_groups);
    SID_log("n_particles=%lld",SID_LOG_COMMENT,n_particles);
    rewind(fp);
    
    // Fill arrays
    n_particles_group=(int     *)SID_malloc(n_groups*sizeof(int));
    group_offsets    =(int     *)SID_malloc(n_groups*sizeof(int));
    particle_ids     =(big_int *)SID_malloc(n_particles*sizeof(big_int));
    i_particle       =0;
    i_group          =0;
    if(n_groups>0)
      group_offsets[0]=0;
    fscanf(fp,"%lld",&n_groups_file);
    while(!feof(fp)){
      for(j_group=0;j_group<n_groups_file;j_group++){
        fscanf(fp,"%d",&(n_particles_group[i_group]));
        for(j_particle=0;j_particle<n_particles_group[i_group];j_particle++,i_particle++)
          fscanf(fp,"%lld",&(particle_ids[i_particle]));
        i_group++;
        if(i_group<n_groups)
          group_offsets[i_group]=group_offsets[i_group-1]+n_particles_group[i_group-1];
      }
      fscanf(fp,"%lld",&n_groups_file);
    }
    fclose(fp);
        
    // Store results
    ADaPS_store(&(plist->data),(void *)(&n_groups),        "n_groups_%s",             ADaPS_SCALAR_INT,   catalog_name);
    ADaPS_store(&(plist->data),(void *)(&n_particles),     "n_particles_%s",          ADaPS_SCALAR_SIZE_T,catalog_name);
    ADaPS_store(&(plist->data),(void *)(n_particles_group),"n_particles_group_%s",    ADaPS_DEFAULT,      catalog_name);
    ADaPS_store(&(plist->data),(void *)(group_offsets),    "particle_offset_group_%s",ADaPS_DEFAULT,      catalog_name);
    ADaPS_store(&(plist->data),(void *)(particle_ids),     "particle_ids_%s",         ADaPS_DEFAULT,      catalog_name);

    SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not open {%s}!",ERROR_IO_READ,filename_ids);

}
