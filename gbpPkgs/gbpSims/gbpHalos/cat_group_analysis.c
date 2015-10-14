#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <sys/stat.h>

int main(int argc, char *argv[]){
  plist_info  plist;
  char        filename_groups_root[256];
  char        filename_snapshot_root[256];
  char        filename_snapshot[256];
  char        filename_number[256];
  char        filename_output_properties[256];
  char        filename_output_properties_dir[256];
  char        filename_output_profiles[256];
  char        filename_output_profiles_dir[256];
  char        filename_output_properties_temp[256];
  char        filename_output_profiles_temp[256];
  char        filename_input_properties[256];
  char        filename_input_properties_dir[256];
  char        filename_input_profiles[256];
  char        filename_input_profiles_dir[256];
  char        filename_input_properties_temp[256];
  char        filename_input_profiles_temp[256];
  char        group_text_prefix[4];
  int         n_groups_process;
  int         n_groups;
  int         n_groups_all;
  int         j_file;
  int         j_file_in;
  int         i_group;
  int         i_file_lo;
  int         i_file_hi;
  int         i_file;
  int         i_file_skip;
  int         i_particle;
  int         j_particle;
  int         i_process;
  int         n_particles;
  int         n_particles_max;
  GBPREAL       *x_array;
  GBPREAL       *y_array;
  GBPREAL       *z_array;
  GBPREAL       *vx_array;
  GBPREAL       *vy_array;
  GBPREAL       *vz_array;
  int        *n_particles_groups_process;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  int        *group_offset;
  size_t      n_particles_in_groups;
  size_t     *ids_snapshot;
  size_t     *ids_groups;
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
  
  FILE       *fp_properties;
  FILE       *fp_profiles;
  FILE       *fp_properties_temp;
  FILE       *fp_profiles_temp;
  cosmo_info *cosmo;
  halo_properties_info  properties;
  halo_profile_info     profile;
  int                   n_files,n_files_i;
  int                   n_files_temp=1;
  int                   flag_process_profiles;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  strcpy(filename_groups_root,argv[1]);
  i_file_lo  =atoi(argv[2]);
  i_file_hi  =atoi(argv[3]);
  i_file_skip=atoi(argv[4]);

  flag_process_profiles=FALSE;

  SID_log("Processing group/subgroup statistics for files #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file_lo,i_file_hi);

  for(i_file=i_file_lo;i_file<=i_file_hi;i_file+=i_file_skip){
    sprintf(filename_number,"%03d", i_file);
    SID_log("Processing file #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);

    for(i_process=0;i_process<2;i_process++){
      switch(i_process){
      case 0:
        sprintf(group_text_prefix,"");
        break;
      case 1:
        sprintf(group_text_prefix,"sub");
        break;
      }
      SID_log("Processing %sgroup files...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
      char filename_base[MAX_FILENAME_LENGTH];
      sprintf(filename_base,"%s_%03d",filename_groups_root,i_file);
      strip_path(filename_base);
      sprintf(filename_input_properties,     "%s_%s.catalog_%sgroups_properties",    filename_groups_root,filename_number,group_text_prefix);
      sprintf(filename_output_properties,    "%s_cat_%s.catalog_%sgroups_properties",filename_groups_root,filename_number,group_text_prefix);
      sprintf(filename_input_properties_dir, "%s",filename_input_properties);
      sprintf(filename_output_properties_dir,"%s",filename_output_properties);
      if(flag_process_profiles){
        sprintf(filename_input_profiles,     "%s_%s.catalog_%sgroups_profiles",    filename_groups_root,filename_number,group_text_prefix);
        sprintf(filename_output_profiles,    "%s_cat_%s.catalog_%sgroups_profiles",filename_groups_root,filename_number,group_text_prefix);
        sprintf(filename_input_profiles_dir, "%s",filename_input_profiles);
        sprintf(filename_output_profiles_dir,"%s",filename_output_profiles);
      }

      // Get the number of files involved
      j_file=0;
      sprintf(filename_input_properties_temp,"%s/%s.catalog_%sgroups_properties.%d",filename_input_properties_dir,filename_base,group_text_prefix,j_file);
      if((fp_properties_temp=fopen(filename_input_properties_temp,"r"))==NULL)
         SID_trap_error("Can not open {%s}",ERROR_IO_OPEN,filename_output_properties_temp);
      fread(&j_file,      sizeof(int),1,fp_properties_temp);
      fread(&n_files,     sizeof(int),1,fp_properties_temp);
      fclose(fp_properties_temp);

      // Concatinate the temporary files together
      fp_properties=fopen(filename_output_properties,"w");
      if(flag_process_profiles)
        fp_profiles  =fopen(filename_output_profiles,  "w");
      for(j_file=0;j_file<n_files;j_file++){

        SID_log("Processing file No. %d of %d...",SID_LOG_OPEN,j_file+1,n_files);
        // Set filenames
        sprintf(filename_input_properties_temp,"%s/%s.catalog_%sgroups_properties.%d",filename_input_properties_dir,filename_base,group_text_prefix,j_file);
        sprintf(filename_input_profiles_temp,  "%s/%s.catalog_%sgroups_profiles.%d",  filename_input_properties_dir,filename_base,group_text_prefix,j_file);
        // Cat properties
        if((fp_properties_temp=fopen(filename_input_properties_temp,"r"))==NULL)
           SID_trap_error("Can't open file {%s}",ERROR_IO_OPEN,filename_output_properties_temp);
        fread(&j_file_in,   sizeof(int),1,fp_properties_temp);
        fread(&n_files_i,   sizeof(int),1,fp_properties_temp);
        fread(&n_groups,    sizeof(int),1,fp_properties_temp);
        fread(&n_groups_all,sizeof(int),1,fp_properties_temp);
        SID_log("Processing properties (%d of %d %sgroups)...",SID_LOG_OPEN,n_groups,n_groups_all,group_text_prefix);
        if(j_file==0){
          fwrite(&j_file,      sizeof(int),1,fp_properties);
          fwrite(&n_files_temp,sizeof(int),1,fp_properties);
          fwrite(&n_groups_all,sizeof(int),1,fp_properties);
          fwrite(&n_groups_all,sizeof(int),1,fp_properties);
        }
        for(i_group=0;i_group<n_groups;i_group++){
          fread( &properties,sizeof(halo_properties_info),1,fp_properties_temp);
          fwrite(&properties,sizeof(halo_properties_info),1,fp_properties);                  
        }
        fclose(fp_properties_temp);
        SID_log("Done.",SID_LOG_CLOSE);

        // Cat profiles
        if(flag_process_profiles){
          fp_profiles_temp=fopen(filename_output_profiles_temp,"r");
          fread(&j_file_in,   sizeof(int),1,fp_profiles_temp);
          fread(&n_files_i,   sizeof(int),1,fp_profiles_temp);
          fread(&n_groups,    sizeof(int),1,fp_profiles_temp);
          fread(&n_groups_all,sizeof(int),1,fp_profiles_temp);
          SID_log("Processing profiles (%d of %d %sgroups)...",SID_LOG_OPEN,n_groups,n_groups_all,group_text_prefix);
          if(j_file==0){
            fwrite(&j_file,      sizeof(int),1,fp_profiles);
            fwrite(&n_files_temp,sizeof(int),1,fp_profiles);
            fwrite(&n_groups_all,sizeof(int),1,fp_profiles);
            fwrite(&n_groups_all,sizeof(int),1,fp_profiles);
          }
          for(i_group=0;i_group<n_groups;i_group++){
            fread( &(profile.n_bins),sizeof(int),1,fp_profiles_temp);
            fwrite(&(profile.n_bins),sizeof(int),1,fp_profiles);
            if(profile.n_bins>0){
              fread( profile.bins,sizeof(halo_profile_bin_info),profile.n_bins,fp_profiles_temp);
              fwrite(profile.bins,sizeof(halo_profile_bin_info),profile.n_bins,fp_profiles);
            }
          }
          fclose(fp_profiles_temp);
          SID_log("Done.",SID_LOG_CLOSE);
        }
        SID_log("Done. (file %3d of %3d)",SID_LOG_CLOSE,j_file,n_files);
      }
      fclose(fp_properties);
      if(flag_process_profiles)
        fclose(fp_profiles);
      SID_log("Done.",SID_LOG_CLOSE);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

