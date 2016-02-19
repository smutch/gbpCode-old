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

int main(int argc, char *argv[]){
  plist_info  plist;
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
  GBPREAL    *x_array;
  GBPREAL    *y_array;
  GBPREAL    *z_array;
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
  char                 *filename_number;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_in_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  int start_snap;
  int stop_snap;
  strcpy(filename_in_root, argv[1]);
  strcpy(filename_out_root,argv[2]);
  start_snap=atoi(argv[3]);
  stop_snap =atoi(argv[4]);

  int offset_size=sizeof(unsigned int);

  SID_log("Converting from root {%s} to root {%s}",SID_LOG_OPEN,filename_in_root,filename_out_root);
  int i_snap;
  for(i_snap=start_snap;i_snap<=stop_snap;i_snap++){
     SID_log("Processing snapshot No. %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap);
     int i_run;
     for (i_run=0;i_run<2;i_run++){
        char filename_in[MAX_FILENAME_LENGTH];
        char filename_out[MAX_FILENAME_LENGTH];
        char group_prefix_text[8];
        switch(i_run){
           case 0:
              sprintf(filename_in, "%s_%03d.catalog_subgroups",filename_in_root, i_snap);
              sprintf(filename_out,"%s_%03d.catalog_subgroups",filename_out_root,i_snap);
              sprintf(group_prefix_text,"sub");
              break;
           case 1:
              sprintf(filename_in, "%s_%03d.catalog_groups",filename_in_root, i_snap);
              sprintf(filename_out,"%s_%03d.catalog_groups",filename_out_root,i_snap);
              sprintf(group_prefix_text,"");
              break;
        }
        FILE *fp_in;
        FILE *fp_out;
        if((fp_in=fopen(filename_in,"r"))==NULL)
           SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_in);
        fp_out=fopen(filename_out,"w");

        // Process n_groups
        int n_groups;
        fread_verify(&n_groups,    sizeof(int),1,fp_in);
        fwrite(&n_groups,   sizeof(int),1,fp_out);
        fwrite(&offset_size,sizeof(int),1,fp_out);

        SID_log("Processing %d %sgroups...",SID_LOG_OPEN,n_groups,group_prefix_text);

        // Process group sizes
        int size;
        for(i_group=0;i_group<n_groups;i_group++){
           fread_verify(&size, sizeof(int),1,fp_in);
           fwrite(&size,sizeof(int),1,fp_out);
        }

        // Process group offsets
        int64_t offset_ll;
        for(i_group=0;i_group<n_groups;i_group++){
           unsigned int offset_in;
           fread_verify(&offset_in, sizeof(unsigned int),1,fp_in);
           if(offset_size==sizeof(unsigned int)){
              fwrite(&offset_in,offset_size,1,fp_out);
           }
           else if(offset_size==sizeof(int64_t)){
              offset_ll=(int64_t)offset_in;
              fwrite(&offset_ll,offset_size,1,fp_out);
           }
           else
              SID_trap_error("Invalid offset size (%d)",ERROR_LOGIC);
        }

        // Process n_subgroups
        if(i_run==1){
           int i_group;
           int n_sub;
           int n_sub_total;
           for(i_group=0,n_sub_total=0;i_group<n_groups;i_group++){
              fread_verify(&n_sub, sizeof(int),1,fp_in);
              fwrite(&n_sub,sizeof(int),1,fp_out);
              n_sub_total+=n_sub;
           }
           SID_log("%d subgroups counted...",SID_LOG_CONTINUE,n_sub_total);
        }
        fclose(fp_in);
        fclose(fp_out);
        SID_log("Done.",SID_LOG_CLOSE);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}

