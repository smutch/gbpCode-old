#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpHalos.h>

#define N_BITS_MIN 1

int main(int argc, char *argv[]){
  
  SID_init(&argc,&argv,NULL,NULL);
  SID_profile_start("make_group_PHKs",SID_PROFILE_NOTMPIENABLED);

  // Fetch user inputs
  char   filename_halos_root[256];
  char   filename_catalog_root[256];
  char   filename_PHKs_root[256];
  double box_size;
  double dx;
  int    i_file_lo_in;
  int    i_file_hi_in;
  int    i_file_skip;
  strcpy(filename_halos_root,  argv[1]);
  strcpy(filename_catalog_root,argv[2]);
  strcpy(filename_PHKs_root,   argv[3]);
  box_size    =atof(argv[4]);
  dx          =atof(argv[5]);
  i_file_lo_in=atoi(argv[6]);
  i_file_hi_in=atoi(argv[7]);
  i_file_skip =atoi(argv[8]);

  int i_file_lo;
  int i_file_hi;
  if(i_file_lo_in<i_file_hi_in){
     i_file_lo=i_file_lo_in;
     i_file_hi=i_file_hi_in;
  }
  else{
     i_file_lo=i_file_hi_in;
     i_file_hi=i_file_lo_in;
  }

  SID_log("Generating group PH keys for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file_lo,i_file_hi);
  for(int i_file=i_file_lo;i_file<=i_file_hi;i_file+=i_file_skip){
    SID_log("Processing file #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);

    // Read group info from the halo catalogs
    plist_info  plist;
    int        *PHK_group      =NULL;
    size_t     *PHK_group_index=NULL;
    char       *filename_number=(char *)SID_malloc(sizeof(char)*10);
    init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
    sprintf(filename_number,"%03d", i_file);
    ADaPS_store(&(plist.data),(void *)filename_number,"read_catalog",ADaPS_DEFAULT);
    read_groups(filename_halos_root,i_file,READ_GROUPS_ALL|READ_GROUPS_MBP_IDS_ONLY,&plist,filename_number);
    int n_groups_all = ((int *)ADaPS_fetch(plist.data,"n_groups_all_%s",filename_number))[0];
    int n_groups     = ((int *)ADaPS_fetch(plist.data,"n_groups_%s",    filename_number))[0];

    // If there's any groups to analyze ...
    int   *n_particles_groups    =NULL;
    size_t n_particles_cumulative=0;
    int    n_bits=0; // Default value if there are no groups
    if(n_groups>0){
       // Fetch the halo sizes
       n_particles_groups=(int *)ADaPS_fetch(plist.data,"n_particles_group_%s",filename_number);

       // Read MBP data from halo catalogs
       SID_log("Reading most-bound-particle positions...",SID_LOG_OPEN);
       halo_properties_info group_properties;
       fp_catalog_info      fp_group_properties;
       double *x_array=(double *)SID_malloc(sizeof(double)*n_groups);
       double *y_array=(double *)SID_malloc(sizeof(double)*n_groups);
       double *z_array=(double *)SID_malloc(sizeof(double)*n_groups);
       fopen_catalog(filename_catalog_root,
                     i_file,
                     READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                     &fp_group_properties);
       if(fp_group_properties.n_halos_total!=n_groups)
         SID_trap_error("Halo counts in group files and catalogs don't match (ie. %d!=%d)",ERROR_LOGIC,fp_group_properties.n_halos_total,n_groups);
       for(int i_group=0;i_group<n_groups;i_group++){
          fread_catalog_file(&fp_group_properties,NULL,&group_properties,NULL,i_group);
          x_array[i_group]=group_properties.position_MBP[0];
          y_array[i_group]=group_properties.position_MBP[1];
          z_array[i_group]=group_properties.position_MBP[2];
          // Enforce periodic BCs
          if(x_array[i_group]<0.)        x_array[i_group]+=box_size;
          if(x_array[i_group]>=box_size) x_array[i_group]-=box_size;
          if(y_array[i_group]<0.)        y_array[i_group]+=box_size;
          if(y_array[i_group]>=box_size) y_array[i_group]-=box_size;
          if(z_array[i_group]<0.)        z_array[i_group]+=box_size;
          if(z_array[i_group]>=box_size) z_array[i_group]-=box_size;
       }
       fclose_catalog(&fp_group_properties);
       SID_log("Done.",SID_LOG_CLOSE);

       // Determine the number of bits to use for the PHKs
       for(n_bits=N_BITS_MIN;(box_size/pow(2.,(double)(n_bits+1)))>dx && n_bits<=20;) n_bits++;
 
       // Compute PHKs
       SID_log("Computing PHKs (using %d bits per dimension)...",SID_LOG_OPEN,n_bits);
       PHK_group=(int *)SID_malloc(sizeof(int)*n_groups);
       for(int i_group=0;i_group<n_groups;i_group++){
         // Compute the key for this group
         PHK_group[i_group]=compute_PHK_from_Cartesian(n_bits,3,(double)x_array[i_group]/box_size,
                                                                (double)y_array[i_group]/box_size,
                                                                (double)z_array[i_group]/box_size);
       }
       SID_free(SID_FARG x_array);
       SID_free(SID_FARG y_array);
       SID_free(SID_FARG z_array);
       SID_log("Done.",SID_LOG_CLOSE);

       // Sort PHKs
       SID_log("Sorting PHKs...",SID_LOG_OPEN);
       merge_sort((void *)PHK_group,n_groups,&PHK_group_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
       SID_log("Done.",SID_LOG_CLOSE);


       // Count the number of particles
       for(int i_group=0;i_group<n_groups;i_group++)
          n_particles_cumulative+=n_particles_groups[PHK_group_index[i_group]];
    }

    // Write results
    SID_log("Writing results for %d groups...",SID_LOG_OPEN,n_groups);
    char filename_output_properties[256];
    sprintf(filename_output_properties,"%s_%s.catalog_PHKs",filename_PHKs_root,filename_number);
    FILE *fp_PHKs=fopen(filename_output_properties,"w");
    fwrite(&n_groups,              sizeof(int),   1,fp_PHKs);
    fwrite(&n_bits,                sizeof(int),   1,fp_PHKs);
    fwrite(&n_particles_cumulative,sizeof(size_t),1,fp_PHKs);
    n_particles_cumulative=0;
    for(int i_group=0;i_group<n_groups;i_group++){
       int index_temp         =(int)PHK_group_index[i_group];
       n_particles_cumulative+=n_particles_groups[index_temp];
       fwrite(&(PHK_group[index_temp]),sizeof(int),   1,fp_PHKs);
       fwrite(&index_temp,             sizeof(int),   1,fp_PHKs);
       fwrite(&n_particles_cumulative, sizeof(size_t),1,fp_PHKs);
    }
    fclose(fp_PHKs);
    SID_log("Done.",SID_LOG_CLOSE);

    // Clean-up
    free_plist(&plist);
    if(n_groups>0){
       SID_free(SID_FARG PHK_group);
       SID_free(SID_FARG PHK_group_index);
    }

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

