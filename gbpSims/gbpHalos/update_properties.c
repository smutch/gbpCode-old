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
  char filename_a_list[MAX_FILENAME_LENGTH];
  char filename_properties_in[MAX_FILENAME_LENGTH];
  char filename_properties_out[MAX_FILENAME_LENGTH];
  char filename_profiles_in[MAX_FILENAME_LENGTH];
  char filename_profiles_out[MAX_FILENAME_LENGTH];
  int start_snap;
  int stop_snap;
  strcpy(filename_in_root, argv[1]);
  strcpy(filename_out_root,argv[2]);
  strcpy(filename_a_list,  argv[3]);
  box_size   =(double)atof(argv[4]);
  start_snap =        atoi(argv[5]);
  stop_snap  =        atoi(argv[6]);

  int offset_size=sizeof(unsigned int);

  SID_log("Converting from root {%s} to root {%s}",SID_LOG_OPEN,filename_in_root,filename_out_root);

  // Read a-list
  double *a_list;
  int     n_a_list;
  FILE   *fp_a_list;
  size_t  line_length=0;
  char   *line=NULL;
  fp_a_list=fopen(filename_a_list,"r");
  n_a_list =count_lines_data(fp_a_list);
  a_list   =(double *)SID_malloc(sizeof(double)*n_a_list);
  for(int i_a_list=0;i_a_list<n_a_list;i_a_list++){
     grab_next_line_data(fp_a_list,&line,&line_length);
     grab_double(line,1,&(a_list[i_a_list]));
  }
  SID_free(SID_FARG line);
  for(int i_a_list=0;i_a_list<n_a_list;i_a_list++)
     fprintf(stderr,"%lf\n",a_list[i_a_list]);
  fclose(fp_a_list);

  int i_snap;
  for(i_snap=start_snap;i_snap<=stop_snap;i_snap++){
     SID_log("Processing snapshot No. %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap);
     double expansion_factor;
     if(i_snap<n_a_list)
        expansion_factor=a_list[i_snap];
     else
        SID_trap_error("Not enough entries in the expansion factor list (ie %d>=%d).",ERROR_LOGIC,i_snap,n_a_list);

     int i_run;
     for (i_run=0;i_run<2;i_run++){
        char filename_in[MAX_FILENAME_LENGTH];
        char filename_out[MAX_FILENAME_LENGTH];
        char group_prefix_text[8];
        switch(i_run){
           case 0:
              sprintf(group_prefix_text,"sub");
              break;
           case 1:
              sprintf(group_prefix_text,"");
              break;
        }
        SID_log("Processing %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,group_prefix_text);
        sprintf(filename_properties_in, "%s_%03d.catalog_%sgroups_properties",filename_in_root, i_snap,group_prefix_text);
        sprintf(filename_properties_out,"%s_%03d.catalog_%sgroups_properties",filename_out_root,i_snap,group_prefix_text);
        sprintf(filename_profiles_in,   "%s_%03d.catalog_%sgroups_profiles",  filename_in_root, i_snap,group_prefix_text);
        sprintf(filename_profiles_out,  "%s_%03d.catalog_%sgroups_profiles",  filename_out_root,i_snap,group_prefix_text);

        // Create filename bases for each dataset
        char filename_properties_in_base[MAX_FILENAME_LENGTH];
        char filename_properties_out_base[MAX_FILENAME_LENGTH];
        char filename_profiles_in_base[MAX_FILENAME_LENGTH];
        char filename_profiles_out_base[MAX_FILENAME_LENGTH];
        strcpy(filename_properties_in_base, filename_properties_in);
        strcpy(filename_properties_out_base,filename_properties_out);
        strcpy(filename_profiles_in_base,   filename_profiles_in);
        strcpy(filename_profiles_out_base,  filename_profiles_out);
        strip_path(filename_properties_in_base);
        strip_path(filename_properties_out_base);
        strip_path(filename_profiles_in_base);
        strip_path(filename_profiles_out_base);

        // Figure out if the properties file(s) are multi-file format
        int   i_file;
        int   n_files_props;
        int   n_props;
        int   n_props_all;
        int   flag_multifile_properties=FALSE;
        FILE *fp_test;
        char  filename_test[MAX_FILENAME_LENGTH];
        sprintf(filename_test,"%s/%s.0",filename_properties_in,filename_properties_in_base);
        if((fp_test=fopen(filename_test,"r"))!=NULL){
           fread(&i_file,       sizeof(int),1,fp_test);
           fread(&n_files_props,sizeof(int),1,fp_test);
           fread(&n_props,      sizeof(int),1,fp_test);
           fread(&n_props_all,  sizeof(int),1,fp_test);
           fclose(fp_test);
           flag_multifile_properties=TRUE;
        }
        else{
           sprintf(filename_test,"%s",filename_properties_in);
           if((fp_test=fopen(filename_test,"r"))!=NULL){
              fread(&i_file,       sizeof(int),1,fp_test);
              fread(&n_files_props,sizeof(int),1,fp_test);
              fread(&n_props,      sizeof(int),1,fp_test);
              fread(&n_props_all,  sizeof(int),1,fp_test);
              fclose(fp_test);
              flag_multifile_properties=FALSE;
              if(n_props!=n_props_all)
                 SID_trap_error("Halo counts don't agree (ie. %d!=%d) in properties file.",ERROR_LOGIC,n_props,n_props_all);
              if(n_files_props!=1)
                 SID_trap_error("Invalid file count (%d) in non-multifile properties file {%s}.",ERROR_LOGIC,n_files_props,filename_test);
           }
           else
              SID_trap_error("Could not open properties dataset.",ERROR_IO_OPEN);
        }
        if(i_file!=0)
           SID_trap_error("Invalid starting file index (%d) in properties file.",ERROR_LOGIC,i_file);

        // Figure out if the profiles file(s) are multi-file format
        int   n_files_profs;
        int   n_profs;
        int   n_profs_all;
        int   flag_multifile_profiles=FALSE;
        sprintf(filename_test,"%s/%s.0",filename_profiles_in,filename_profiles_in_base);
        if((fp_test=fopen(filename_test,"r"))!=NULL){
           int i_file;
           fread(&i_file,       sizeof(int),1,fp_test);
           fread(&n_files_profs,sizeof(int),1,fp_test);
           fread(&n_profs,      sizeof(int),1,fp_test);
           fread(&n_profs_all,  sizeof(int),1,fp_test);
           fclose(fp_test);
           flag_multifile_profiles=TRUE;
        }
        else{
           sprintf(filename_test,"%s",filename_profiles_in);
           if((fp_test=fopen(filename_test,"r"))!=NULL){
              fread(&i_file,       sizeof(int),1,fp_test);
              fread(&n_files_profs,sizeof(int),1,fp_test);
              fread(&n_profs,      sizeof(int),1,fp_test);
              fread(&n_profs_all,  sizeof(int),1,fp_test);
              fclose(fp_test);
              flag_multifile_profiles=FALSE;
              if(n_profs!=n_profs_all)
                 SID_trap_error("Halo counts don't agree (ie. %d!=%d) in profiles file.",ERROR_LOGIC,n_profs,n_profs_all);
              if(n_files_profs!=1)
                 SID_trap_error("Invalid file count (%d) in non-multifile profiles file.",ERROR_LOGIC,n_files_profs);
           }  
           else
              SID_trap_error("Could not open profiles dataset.",ERROR_IO_OPEN);
        }
        if(i_file!=0)
           SID_trap_error("Invalid starting file index (%d) in profiles file.",ERROR_LOGIC,i_file);

        // Check that the halo counts in the properties and profiles datasets agree
        if(n_profs_all!=n_props_all)
           SID_trap_error("The properties and profiles halo counts don't agree (ie. %d!=%d)",ERROR_LOGIC,n_profs_all,n_props_all);

        // Loop over all the halos
        FILE *fp_properties_read =NULL;
        FILE *fp_properties_write=NULL;
        FILE *fp_profiles_read   =NULL;
        FILE *fp_profiles_write  =NULL;
        int   i_halo_props_last;
        int   i_halo_profs_last;
        int   i_halo;
        int   i_file_props=0;
        int   i_file_profs=0;
        if(n_props_all==0){
           char filename_read[MAX_FILENAME_LENGTH];
           char filename_write[MAX_FILENAME_LENGTH];
           // Copy empty property files
           if(i_file_props==0 && flag_multifile_properties)
              mkdir(filename_properties_out,02755);
           for(i_file_props=0;i_file_props<n_files_props;i_file_props++){
              // Open files
              if(flag_multifile_properties){
                 sprintf(filename_read, "%s/%s.%d",filename_properties_in, filename_properties_in_base, i_file_props);
                 sprintf(filename_write,"%s/%s.%d",filename_properties_out,filename_properties_out_base,i_file_props);
              }
              else{
                 if(i_file_props>0)
                    SID_trap_error("Error looping over non-multifile.",ERROR_LOGIC);
                 sprintf(filename_read, "%s",filename_properties_in);
                 sprintf(filename_write,"%s",filename_properties_out);
              }
              if((fp_properties_read =fopen(filename_read, "r"))==NULL)
                 SID_trap_error("Could not open file {%s} for reading.",ERROR_IO_OPEN,filename_read);
              if((fp_properties_write=fopen(filename_write,"w"))==NULL)
                 SID_trap_error("Could not open file {%s} for writing.",ERROR_IO_OPEN,filename_write);
              // Read header
              fread(&i_file,       sizeof(int),1,fp_properties_read);
              fread(&n_files_props,sizeof(int),1,fp_properties_read);
              fread(&n_props,      sizeof(int),1,fp_properties_read);
              fread(&n_props_all,  sizeof(int),1,fp_properties_read);
              // Write header
              fwrite(&i_file,       sizeof(int),1,fp_properties_write);
              fwrite(&n_files_props,sizeof(int),1,fp_properties_write);
              fwrite(&n_props,      sizeof(int),1,fp_properties_write);
              fwrite(&n_props_all,  sizeof(int),1,fp_properties_write);
              // Close files
              fclose(fp_properties_read);
              fclose(fp_properties_write);
           }

           // Copy empty profile files
           if(i_file_profs==0 && flag_multifile_profiles)
              mkdir(filename_profiles_out,02755);
           for(i_file_profs=0;i_file_profs<n_files_profs;i_file_profs++){
              // Open files
              if(flag_multifile_profiles){
                 sprintf(filename_read, "%s/%s.%d",filename_profiles_in, filename_profiles_in_base, i_file_profs);
                 sprintf(filename_write,"%s/%s.%d",filename_profiles_out,filename_profiles_out_base,i_file_profs);
              }
              else{
                 if(i_file_profs>0)
                    SID_trap_error("Error looping over non-multifile.",ERROR_LOGIC);
                 sprintf(filename_read, "%s",filename_profiles_in);
                 sprintf(filename_write,"%s",filename_profiles_out);
              }
              if((fp_profiles_read =fopen(filename_read, "r"))==NULL)
                 SID_trap_error("Could not open file {%s} for reading.",ERROR_IO_OPEN,filename_read);
              if((fp_profiles_write=fopen(filename_write,"w"))==NULL)
                 SID_trap_error("Could not open file {%s} for writing.",ERROR_IO_OPEN,filename_write);
              // Read header
              fread(&i_file,       sizeof(int),1,fp_profiles_read);
              fread(&n_files_profs,sizeof(int),1,fp_profiles_read);
              fread(&n_profs,      sizeof(int),1,fp_profiles_read);
              fread(&n_profs_all,  sizeof(int),1,fp_profiles_read);
              // Write header
              fwrite(&i_file,       sizeof(int),1,fp_profiles_write);
              fwrite(&n_files_profs,sizeof(int),1,fp_profiles_write);
              fwrite(&n_profs,      sizeof(int),1,fp_profiles_write);
              fwrite(&n_profs_all,  sizeof(int),1,fp_profiles_write);
              // Close files
              fclose(fp_profiles_read);
              fclose(fp_profiles_write);
           }  
        }
        else{
           for(i_halo=0,i_halo_props_last=0,i_halo_profs_last=0,i_file_props=0,i_file_profs=0;i_halo<n_props_all;i_halo++){
              // Open properties files and read/write their headers
              if(i_halo>=i_halo_props_last){
                 // Open files
                 char filename_read[MAX_FILENAME_LENGTH];
                 char filename_write[MAX_FILENAME_LENGTH];
                 if(flag_multifile_properties){
                    if(i_file_props==0)
                       mkdir(filename_properties_out,02755);
                    sprintf(filename_read, "%s/%s.%d",filename_properties_in, filename_properties_in_base, i_file_props);
                    sprintf(filename_write,"%s/%s.%d",filename_properties_out,filename_properties_out_base,i_file_props);
                    i_file_props++;
                 }
                 else{
                    sprintf(filename_read, "%s",filename_properties_in);
                    sprintf(filename_write,"%s",filename_properties_out);
                 }
                 if(fp_properties_read!=NULL)
                    fclose(fp_properties_read);
                 if(fp_properties_write!=NULL)
                    fclose(fp_properties_write);
                 if((fp_properties_read =fopen(filename_read, "r"))==NULL)
                    SID_trap_error("Could not open file {%s} for reading.",ERROR_IO_OPEN,filename_read);
                 if((fp_properties_write=fopen(filename_write,"w"))==NULL)
                    SID_trap_error("Could not open file {%s} for writing.",ERROR_IO_OPEN,filename_write);

                 // Read header
                 fread(&i_file,       sizeof(int),1,fp_properties_read);
                 fread(&n_files_props,sizeof(int),1,fp_properties_read);
                 fread(&n_props,      sizeof(int),1,fp_properties_read);
                 fread(&n_props_all,  sizeof(int),1,fp_properties_read);

                 // Write header
                 fwrite(&i_file,       sizeof(int),1,fp_properties_write);
                 fwrite(&n_files_props,sizeof(int),1,fp_properties_write);
                 fwrite(&n_props,      sizeof(int),1,fp_properties_write);
                 fwrite(&n_props_all,  sizeof(int),1,fp_properties_write);

                 // Increment counter
                 i_halo_props_last+=n_props;
              }

              // Open profiles files and read/write their headers
              if(i_halo>=i_halo_profs_last){
                 // Open files
                 char filename_read[MAX_FILENAME_LENGTH];
                 char filename_write[MAX_FILENAME_LENGTH];
                 if(flag_multifile_profiles){
                    if(i_file_profs==0)
                       mkdir(filename_profiles_out,02755);
                    sprintf(filename_read, "%s/%s.%d",filename_profiles_in, filename_profiles_in_base, i_file_profs);
                    sprintf(filename_write,"%s/%s.%d",filename_profiles_out,filename_profiles_out_base,i_file_profs);
                    i_file_profs++;
                 }
                 else{
                    sprintf(filename_read, "%s",filename_profiles_in);
                    sprintf(filename_write,"%s",filename_profiles_out);
                 }
                 if(fp_profiles_read!=NULL)
                    fclose(fp_profiles_read);
                 if(fp_profiles_write!=NULL)
                    fclose(fp_profiles_write);
                 if((fp_profiles_read =fopen(filename_read, "r"))==NULL)
                    SID_trap_error("Could not open file {%s} for reading.",ERROR_IO_OPEN,filename_read);
                 if((fp_profiles_write=fopen(filename_write,"w"))==NULL)
                    SID_trap_error("Could not open file {%s} for writing.",ERROR_IO_OPEN,filename_write);

                 // Read header
                 fread(&i_file,       sizeof(int),1,fp_profiles_read);
                 fread(&n_files_profs,sizeof(int),1,fp_profiles_read);
                 fread(&n_profs,      sizeof(int),1,fp_profiles_read);
                 fread(&n_profs_all,  sizeof(int),1,fp_profiles_read);

                 // Write header
                 fwrite(&i_file,       sizeof(int),1,fp_profiles_write);
                 fwrite(&n_files_profs,sizeof(int),1,fp_profiles_write);
                 fwrite(&n_profs,      sizeof(int),1,fp_profiles_write);
                 fwrite(&n_profs_all,  sizeof(int),1,fp_profiles_write);

                 // Increment counter
                 i_halo_profs_last+=n_profs;
              }

              // Read profiles
              int                   n_bins;
              halo_profile_bin_info bins[MAX_PROFILE_BINS];
              fread(&n_bins,sizeof(int),                  1,     fp_profiles_read);
              fread(&bins,  sizeof(halo_profile_bin_info),n_bins,fp_profiles_read);

              // Modify profiles
              /*
              int i_bin;
              int n_particles_cumulative=0;
              for(i_bin=0;i_bin<n_bins;i_bin++){
                 n_particles_cumulative+=bins[i_bin].n_particles;
                 bins[i_bin].spin[0]   /=(float)n_particles_cumulative;
                 bins[i_bin].spin[1]   /=(float)n_particles_cumulative;
                 bins[i_bin].spin[2]   /=(float)n_particles_cumulative;
              }
              */

              // Write profiles
              fwrite(&n_bins,sizeof(int),                  1,     fp_profiles_write);
              fwrite(&bins,  sizeof(halo_profile_bin_info),n_bins,fp_profiles_write);

              // Read profiles
              halo_properties_info  properties;
              fread(&properties,sizeof(halo_properties_info),1,fp_properties_read);

              // Modify properties
              /*
              const gsl_interp_type *interp_type;
              interp_info *vir_interpolate;
              double       r_interp[MAX_PROFILE_BINS];
              double       y_interp[MAX_PROFILE_BINS];
              if(n_bins>9)
                 interp_type=gsl_interp_cspline;
              else
                 interp_type=gsl_interp_linear;
              interp_type=gsl_interp_linear;
              for(i_bin=0;i_bin<n_bins;i_bin++)
                 r_interp[i_bin]=(double)bins[i_bin].r_max;
              for(i_bin=0;i_bin<n_bins;i_bin++)
                 y_interp[i_bin]=(double)bins[i_bin].spin[0];
              init_interpolate(r_interp,y_interp,n_bins,interp_type,&vir_interpolate);
              properties.spin[0]=(float)interpolate(vir_interpolate,properties.R_vir);
              free_interpolate(SID_FARG vir_interpolate,NULL);
              for(i_bin=0;i_bin<n_bins;i_bin++)
                 y_interp[i_bin]=(double)bins[i_bin].spin[1];
              init_interpolate(r_interp,y_interp,n_bins,interp_type,&vir_interpolate);
              properties.spin[1]=(float)interpolate(vir_interpolate,properties.R_vir);
              free_interpolate(SID_FARG vir_interpolate,NULL);
              for(i_bin=0;i_bin<n_bins;i_bin++)
                 y_interp[i_bin]=(double)bins[i_bin].spin[2];
              init_interpolate(r_interp,y_interp,n_bins,interp_type,&vir_interpolate);
              properties.spin[2]=(float)interpolate(vir_interpolate,properties.R_vir);
              free_interpolate(SID_FARG vir_interpolate,NULL);
              */
              if(properties.R_halo==properties.R_vir){
                 int i_profile=n_bins-1;
                 //  ... COM positions ...
                 properties.position_COM[0]=(double)(bins[i_profile].position_COM[0])/expansion_factor;
                 properties.position_COM[1]=(double)(bins[i_profile].position_COM[1])/expansion_factor;
                 properties.position_COM[2]=(double)(bins[i_profile].position_COM[2])/expansion_factor;
                 //  ... M_vir ...
                 properties.M_vir=(double)bins[i_profile].M_r;
                 //  ... sigma_v ...
                 properties.sigma_v=(float)bins[i_profile].sigma_tot;
                 //   ... spin ...
                 properties.spin[0]=(float)bins[i_profile].spin[0];
                 properties.spin[1]=(float)bins[i_profile].spin[1];
                 properties.spin[2]=(float)bins[i_profile].spin[2];
                 //  ... triaxial axes ratios ...
                 properties.q_triaxial=(float)bins[i_profile].q_triaxial;
                 properties.s_triaxial=(float)bins[i_profile].s_triaxial;
                 // ... shape eigen vectors ...
                 for(int i=0;i<3;i++){
                   double norm;
                   for(int j=0;j<3;j++)
                     properties.shape_eigen_vectors[i][j]=(float)bins[i_profile].shape_eigen_vectors[i][j];
                   norm=sqrt(properties.shape_eigen_vectors[i][0]*properties.shape_eigen_vectors[i][0]+
                             properties.shape_eigen_vectors[i][1]*properties.shape_eigen_vectors[i][1]+
                             properties.shape_eigen_vectors[i][2]*properties.shape_eigen_vectors[i][2]);
                   for(int j=0;j<3;j++)
                     properties.shape_eigen_vectors[i][j]/=norm;
                 }
              }

              // Enforce periodic box on COM position
              properties.position_COM[0]+=properties.position_MBP[0];
              properties.position_COM[1]+=properties.position_MBP[1];
              properties.position_COM[2]+=properties.position_MBP[2];
              if(properties.position_COM[0]< box_size) properties.position_COM[0]+=box_size;
              if(properties.position_COM[1]< box_size) properties.position_COM[1]+=box_size;
              if(properties.position_COM[2]< box_size) properties.position_COM[2]+=box_size;
              if(properties.position_COM[0]>=box_size) properties.position_COM[0]-=box_size;
              if(properties.position_COM[1]>=box_size) properties.position_COM[1]-=box_size;
              if(properties.position_COM[2]>=box_size) properties.position_COM[2]-=box_size;

              // Write properties
              fwrite(&properties,sizeof(halo_properties_info),1,fp_properties_write);
           }

           // Perform final close
           if(fp_properties_read!=NULL)
              fclose(fp_properties_read);
           if(fp_properties_write!=NULL)
              fclose(fp_properties_write);
           if(fp_profiles_read!=NULL)
              fclose(fp_profiles_read);
           if(fp_profiles_write!=NULL)
              fclose(fp_profiles_write);
        }

        SID_log("Done.",SID_LOG_CLOSE);
     }
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_free(SID_FARG a_list);
  SID_exit(ERROR_NONE);
}

