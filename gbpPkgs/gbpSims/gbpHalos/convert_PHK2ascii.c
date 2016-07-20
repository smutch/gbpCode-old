#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpHalos.h>

#define N_BITS_MIN 1

int main(int argc, char *argv[]){
  plist_info  plist;
  char        filename_halos_root[256];
  char        filename_catalog_root[256];
  char        filename_PHKs_root[256];
  char        filename_snapshot_root[256];
  char        filename_snapshot[256];
  char        filename_output_properties_dir[256];
  char        filename_output_properties[256];
  char        filename_output_profiles_dir[256];
  char        filename_output_profiles[256];
  char        filename_output_properties_temp[256];
  char        filename_output_profiles_temp[256];
  char        group_text_prefix[4];
  int         n_groups_process;
  int         n_groups;
  int         n_groups_all;
  int         i_rank;
  int         i_group;
  int         i_file_lo_in;
  int         i_file_lo;
  int         i_file_hi_in;
  int         i_file_hi;
  int         i_file;
  int         i_file_skip;
  int         i_particle;
  int         j_particle;
  int         i_process;
  int         n_particles;
  int         n_particles_max;
  size_t      n_particles_cumulative;
  double     *x_array;
  double     *y_array;
  double     *z_array;
  size_t     *ids_particles;
  size_t     *ids_particles_index;
  size_t     *ids_groups;
  int        *n_particles_groups_process;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  int        *group_offset;
  size_t      n_particles_in_groups;
  size_t     *ids_snapshot;
  size_t     *ids_sort_index;
  size_t     *ids_snapshot_sort_index;
  size_t     *ids_groups_sort_index;
  size_t      n_particles_snapshot;
  int         i_value;
  double      h_Hubble;
  double      Omega_M;
  double      Omega_b;
  double      Omega_Lambda;
  double      f_gas;
  double      Omega_k;
  double      sigma_8;
  double      n_spec;
  double      redshift;
  double      box_size;
  double      dx;
  double      particle_mass;
  
  int         r_val;
  struct      stat file_stats;
  size_t      n_bytes;
  size_t      n_bytes_buffer;
  void       *buffer;
  
  FILE       *fp_PHKs;
  FILE       *fp_profiles;
  FILE       *fp_PHKs_temp;
  FILE       *fp_profiles_temp;
  cosmo_info *cosmo;
  halo_properties_info  properties;
  halo_profile_info     profile;
  int                   n_temp;
  int                   n_truncated;
  int                   largest_truncated;
  int                   largest_truncated_local;
  int                   n_bits=N_BITS_MIN;
  char                 *filename_number;

  SID_init(&argc,&argv,NULL,NULL);
  SID_profile_start("make_group_PHKs",SID_PROFILE_NOTMPIENABLED);

  // Fetch user inputs
  strcpy(filename_PHKs_root,argv[1]);
  i_file_lo_in        =atoi(argv[2]);
  i_file_hi_in        =atoi(argv[3]);
  i_file_skip         =atoi(argv[4]);

  if(i_file_lo_in<i_file_hi_in){
     i_file_lo=i_file_lo_in;
     i_file_hi=i_file_hi_in;
  }
  else{
     i_file_lo=i_file_hi_in;
     i_file_hi=i_file_lo_in;
  }

  SID_log("Converting group PH keys to ascii for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file_lo,i_file_hi);
  for(i_file=i_file_hi;i_file>=i_file_lo;i_file-=i_file_skip){
    SID_log("Processing file #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);

    // Open files
    FILE *fp_in =NULL;
    FILE *fp_out=NULL;
    char  filename_in[MAX_FILENAME_LENGTH];
    char  filename_out[MAX_FILENAME_LENGTH];
    sprintf(filename_in, "%s_%03d.catalog_PHKs",      filename_PHKs_root,i_file);
    sprintf(filename_out,"%s_%03d.catalog_PHKs.ascii",filename_PHKs_root,i_file);
    if((fp_in=fopen(filename_in,"r"))==NULL)
       SID_trap_error("Could not open input file {%s}",ERROR_IO_OPEN,filename_in);
    if((fp_out=fopen(filename_out,"w"))==NULL)
       SID_trap_error("Could not open output file {%s}",ERROR_IO_OPEN,filename_out);

    // Read/write the header
    fread_verify(&n_groups,              sizeof(int),   1,fp_in);
    fread_verify(&n_bits,                sizeof(int),   1,fp_in);
    fread_verify(&n_particles_cumulative,sizeof(size_t),1,fp_in);
    fprintf(fp_out,"# n_groups   = %d\n",n_groups);
    fprintf(fp_out,"# n_bits     = %d\n",n_bits);
    fprintf(fp_out,"# n_particles= %d\n",n_particles_cumulative);

    // Read/write each group in turn
    int    PHK_group;
    int    index_temp;
    size_t n_particles_cumulative;
    for(i_group=0;i_group<n_groups;i_group++){
       fread_verify(&PHK_group,             sizeof(int),   1,fp_in);
       fread_verify(&index_temp,            sizeof(int),   1,fp_in);
       fread_verify(&n_particles_cumulative,sizeof(size_t),1,fp_in);
       fprintf(fp_out,"%6d %6d %6d %lld\n",i_group,PHK_group,index_temp,n_particles_cumulative);
    }
    fclose(fp_in);
    fclose(fp_out);

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

