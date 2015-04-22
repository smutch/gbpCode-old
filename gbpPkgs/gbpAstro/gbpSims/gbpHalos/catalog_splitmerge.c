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
  char        filename_output_SO_dir[256];
  char        filename_output_SO[256];
  char        filename_output_properties_temp[256];
  char        filename_output_profiles_temp[256];
  char        filename_output_SO_temp[256];
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
  char filename_properties_in[MAX_FILENAME_LENGTH];
  char filename_properties_out[MAX_FILENAME_LENGTH];
  char filename_profiles_in[MAX_FILENAME_LENGTH];
  char filename_profiles_out[MAX_FILENAME_LENGTH];
  char filename_SO_in[MAX_FILENAME_LENGTH];
  char filename_SO_out[MAX_FILENAME_LENGTH];
  strcpy(filename_in_root,    argv[1]);
  strcpy(filename_out_root,   argv[2]);
  int start_snap       = atoi(argv[3]);
  int stop_snap        = atoi(argv[4]);
  int n_files_out_in   = atoi(argv[5]);

  // Process each snapshot in turn
  SID_log("Rewriting catalogs from root {%s} to root {%s} with %d files...",SID_LOG_OPEN,filename_in_root,filename_out_root,n_files_out_in);
  for(int i_snap=start_snap;i_snap<=stop_snap;i_snap++){
     SID_log("Processing snapshot No. %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap);

     // Loop once for subgroups and once for groups
     for (int i_run=0;i_run<2;i_run++){
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
        sprintf(filename_SO_in,         "%s_%03d.catalog_%sgroups_SO",        filename_in_root, i_snap,group_prefix_text);
        sprintf(filename_SO_out,        "%s_%03d.catalog_%sgroups_SO",        filename_out_root,i_snap,group_prefix_text);

        // Create filename bases for each dataset
        char filename_properties_in_base[MAX_FILENAME_LENGTH];
        char filename_properties_out_base[MAX_FILENAME_LENGTH];
        char filename_profiles_in_base[MAX_FILENAME_LENGTH];
        char filename_profiles_out_base[MAX_FILENAME_LENGTH];
        char filename_SO_in_base[MAX_FILENAME_LENGTH];
        char filename_SO_out_base[MAX_FILENAME_LENGTH];
        strcpy(filename_properties_in_base, filename_properties_in);
        strcpy(filename_properties_out_base,filename_properties_out);
        strcpy(filename_profiles_in_base,   filename_profiles_in);
        strcpy(filename_profiles_out_base,  filename_profiles_out);
        strcpy(filename_SO_in_base,         filename_SO_in);
        strcpy(filename_SO_out_base,        filename_SO_out);
        strip_path(filename_properties_in_base);
        strip_path(filename_properties_out_base);
        strip_path(filename_profiles_in_base);
        strip_path(filename_profiles_out_base);
        strip_path(filename_SO_in_base);
        strip_path(filename_SO_out_base);

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
           fread_verify(&i_file,       sizeof(int),1,fp_test);
           fread_verify(&n_files_props,sizeof(int),1,fp_test);
           fread_verify(&n_props,      sizeof(int),1,fp_test);
           fread_verify(&n_props_all,  sizeof(int),1,fp_test);
           fclose(fp_test);
           flag_multifile_properties=TRUE;
        }
        else{
           sprintf(filename_test,"%s",filename_properties_in);
           if((fp_test=fopen(filename_test,"r"))!=NULL){
              fread_verify(&i_file,       sizeof(int),1,fp_test);
              fread_verify(&n_files_props,sizeof(int),1,fp_test);
              fread_verify(&n_props,      sizeof(int),1,fp_test);
              fread_verify(&n_props_all,  sizeof(int),1,fp_test);
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
           fread_verify(&i_file,       sizeof(int),1,fp_test);
           fread_verify(&n_files_profs,sizeof(int),1,fp_test);
           fread_verify(&n_profs,      sizeof(int),1,fp_test);
           fread_verify(&n_profs_all,  sizeof(int),1,fp_test);
           fclose(fp_test);
           flag_multifile_profiles=TRUE;
        }
        else{
           sprintf(filename_test,"%s",filename_profiles_in);
           if((fp_test=fopen(filename_test,"r"))!=NULL){
              fread_verify(&i_file,       sizeof(int),1,fp_test);
              fread_verify(&n_files_profs,sizeof(int),1,fp_test);
              fread_verify(&n_profs,      sizeof(int),1,fp_test);
              fread_verify(&n_profs_all,  sizeof(int),1,fp_test);
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
        int n_halos_all=n_props_all;

        // SO files are defined only for groups
        int n_files_SO;
        int n_SO;
        int n_SO_all;
        int flag_multifile_SO=FALSE;
        int flag_SO_present  =TRUE;
        if(i_run==1){
           // Figure out if the SO file(s) are multi-file format
           sprintf(filename_test,"%s/%s.0",filename_SO_in,filename_SO_in_base);
           if((fp_test=fopen(filename_test,"r"))!=NULL){
              int i_file;
              fread_verify(&i_file,    sizeof(int),1,fp_test);
              fread_verify(&n_files_SO,sizeof(int),1,fp_test);
              fread_verify(&n_SO,      sizeof(int),1,fp_test);
              fread_verify(&n_SO_all,  sizeof(int),1,fp_test);
              fclose(fp_test);
              flag_multifile_SO=TRUE;
           }
           else{
              sprintf(filename_test,"%s",filename_SO_in);
              if((fp_test=fopen(filename_test,"r"))!=NULL){
                 fread_verify(&i_file,    sizeof(int),1,fp_test);
                 fread_verify(&n_files_SO,sizeof(int),1,fp_test);
                 fread_verify(&n_SO,      sizeof(int),1,fp_test);
                 fread_verify(&n_SO_all,  sizeof(int),1,fp_test);
                 fclose(fp_test);
                 flag_multifile_SO=FALSE;
                 if(n_SO!=n_SO_all)
                    SID_trap_error("Halo counts don't agree (ie. %d!=%d) in SO file.",ERROR_LOGIC,n_SO,n_SO_all);
                 if(n_files_SO!=1)
                    SID_trap_error("Invalid file count (%d) in non-multifile SO file.",ERROR_LOGIC,n_files_SO);
              }
              // Not all snapshots are garanteed to have an SO dataset.  Continue if one isn't found
              else{
                 flag_SO_present=FALSE;
                 SID_log("Could not open SO dataset.",SID_LOG_COMMENT);
              }
           }
           if(i_file!=0)
              SID_trap_error("Invalid starting file index (%d) in SO file.",ERROR_LOGIC,i_file);

           // Check that the halo counts in the properties and profiles datasets agree
           if(n_SO_all!=n_props_all && flag_SO_present)
              SID_trap_error("The properties and SO halo counts don't agree (ie. %d!=%d)",ERROR_LOGIC,n_SO_all,n_props_all);
        }

        // If there are fewer than 1000 halos, don't bother
        //    writting to a multi-file
        int n_files_out=n_files_out_in;
        if(n_halos_all<MAX(n_files_out_in,1000))
           n_files_out=1;

        // Perform rewrites
        int n_rewrite=2;
        if(i_run==1 && flag_SO_present) n_rewrite=3;
        for(int i_rewrite=2;i_rewrite<n_rewrite;i_rewrite++){
           char *buffer  =NULL;
           FILE *fp_read =NULL;
           FILE *fp_write=NULL;
           int   i_file_read  =0;
           int   i_file_write =0;
           int   i_halo_read  =0; // file index
           int   i_halo_write =0; // file index
           int   j_halo_read  =0; // catalog_index
           int   j_halo_write =0; // catalog index
           int   n_halos_read =0;
           int   n_halos_write=0;
           int   n_files_rewrite;
           int   n_items_all;
           int   flag_multifile;
           char *filename_items_in;
           char *filename_items_in_base;
           char *filename_items_out;
           char *filename_items_out_base;
           if(i_rewrite==0){
              SID_log("Processing properties...",SID_LOG_OPEN|SID_LOG_TIMER);
              buffer                 =(char *)SID_malloc(sizeof(halo_properties_info));
              n_files_rewrite        =n_files_props;
              n_items_all            =n_props_all;
              flag_multifile         =flag_multifile_properties;
              filename_items_in      =filename_properties_in;
              filename_items_in_base =filename_properties_in_base;
              filename_items_out     =filename_properties_out;
              filename_items_out_base=filename_properties_out_base;
           }
           else if(i_rewrite==1){
              SID_log("Processing profiles...",SID_LOG_OPEN|SID_LOG_TIMER);
              buffer                 =(char *)SID_malloc(sizeof(halo_profile_bin_info));
              n_files_rewrite        =n_files_profs;
              n_items_all            =n_profs_all;
              flag_multifile         =flag_multifile_profiles;
              filename_items_in      =filename_profiles_in;
              filename_items_in_base =filename_profiles_in_base;
              filename_items_out     =filename_profiles_out;
              filename_items_out_base=filename_profiles_out_base;
           }
           else if(i_rewrite==2){
              SID_log("Processing SO files...",SID_LOG_OPEN|SID_LOG_TIMER);
              buffer                 =(char *)SID_malloc(6*sizeof(float));
              n_files_rewrite        =n_files_SO;
              n_items_all            =n_SO_all;
              flag_multifile         =flag_multifile_SO;
              filename_items_in      =filename_SO_in;
              filename_items_in_base =filename_SO_in_base;
              filename_items_out     =filename_SO_out;
              filename_items_out_base=filename_SO_out_base;
           }
           for(int i_halo=0;i_halo<n_halos_all;i_halo++){
              // Open a new input file if need-be
              while(i_halo_read>=n_halos_read && i_file_read<n_files_rewrite){
                 if(!flag_multifile && i_file_read>0)
                    SID_trap_error("Trying to open a second file in a non-multifile dataset.",ERROR_LOGIC);
                 if(flag_multifile)
                    sprintf(filename_in,"%s/%s.%d",filename_items_in,filename_items_in_base,i_file_read);
                 else
                    sprintf(filename_in,"%s",filename_items_in);
                 if(fp_read!=NULL)
                    fclose(fp_read);
                 if((fp_read=fopen(filename_in,"r"))!=NULL){
                    int n_files_in;
                    fread_verify(&i_file,      sizeof(int),1,fp_read);
                    fread_verify(&n_files_in,  sizeof(int),1,fp_read);
                    fread_verify(&n_halos_read,sizeof(int),1,fp_read);
                    fread_verify(&n_items_all, sizeof(int),1,fp_read);
                    if(n_files_in!=n_files_rewrite)
                       SID_trap_error("File counts are not consistant (ie. %d!=%d).",ERROR_LOGIC,n_files_in,n_files_rewrite);
                 }
                 else
                    SID_trap_error("Could not open {%s}.",ERROR_IO_OPEN,filename_in);
                 if(i_file!=i_file_read)
                    SID_trap_error("Invalid file index in (ie. %d!=%d).",ERROR_LOGIC,i_file,i_file_read);
                 if(n_items_all!=n_halos_all)
                    SID_trap_error("Invalid total halo count in {%s} (ie. %d!=%d).",ERROR_LOGIC,filename_in,n_items_all,n_halos_all);
                 i_halo_read=0;
                 i_file_read++;
              }
              // Open a new output file if need-be
              if(i_halo_write>=n_halos_write){
                 if(n_files_out==0 && i_file_write>0)
                    SID_trap_error("Trying to create a second file in a non-multifile dataset.",ERROR_LOGIC);
                 if(n_files_out>1){
                    if(i_file_write==0)
                       mkdir(filename_items_out,02755);
                    sprintf(filename_out,"%s/%s.%d",filename_items_out,filename_items_out_base,i_file_write);
                 }
                 else
                    sprintf(filename_out,"%s",filename_items_out);
                 if(i_file_write==(n_files_out-1))
                    n_halos_write=n_halos_all-i_halo_write;
                 else
                    n_halos_write=(int)((float)(i_file_write+1)*(float)n_halos_all/(float)n_files_out)-j_halo_write;
                 int n_halos_left=n_halos_all-j_halo_write;
                 if(n_halos_write>n_halos_left)
                    n_halos_write=n_halos_left;
                 if(fp_write!=NULL)
                    fclose(fp_write);
                 if((fp_write=fopen(filename_out,"w"))!=NULL){
                    fwrite(&i_file_write, sizeof(int),1,fp_write);
                    fwrite(&n_files_out,  sizeof(int),1,fp_write);
                    fwrite(&n_halos_write,sizeof(int),1,fp_write);
                    fwrite(&n_halos_all,  sizeof(int),1,fp_write);
                 }
                 else
                    SID_trap_error("Could not create properties file {%s}.",ERROR_IO_OPEN,filename_out);
                 i_halo_write=0;
                 i_file_write++;
              }
              switch(i_rewrite){
                 // Rewrite properties
                 case 0:
                    fread_verify(buffer,sizeof(halo_properties_info),1,fp_read);
                    fwrite      (buffer,sizeof(halo_properties_info),1,fp_write);
                    break;
                 // Rewrite profiles
                 case 1:
                    int n_bins;
                    fread_verify(&n_bins,sizeof(int),1,fp_read);
                    fwrite      (&n_bins,sizeof(int),1,fp_write);
                    for(int i_bin=0;i_bin<n_bins;i_bin++){
                       fread_verify(buffer,sizeof(halo_profile_bin_info),1,fp_read);
                       fwrite      (buffer,sizeof(halo_profile_bin_info),1,fp_write);
                    }
                    break;
                 // Rewrite SOs
                 case 2:
                    fread_verify(buffer,sizeof(float),6,fp_read);
                    fwrite      (buffer,sizeof(float),6,fp_write);
                    break;
              }
              i_halo_read++;
              j_halo_read++;
              i_halo_write++;
              j_halo_write++;
           } // i_halo
           SID_free(SID_FARG buffer);
           if(fp_read!=NULL)
              fclose(fp_read);
           if(fp_write!=NULL)
              fclose(fp_write);
           fp_read =NULL;
           fp_write=NULL;
           // Sanity check
           if(j_halo_read!=n_halos_all)
              SID_trap_error("The proper number of halos was not read (ie. %d!=%d).",ERROR_IO_OPEN,j_halo_read,n_halos_all);
           if(j_halo_write!=n_halos_all)
              SID_trap_error("The proper number of halos was not written (ie. %d!=%d).",ERROR_IO_OPEN,j_halo_write,n_halos_all);
           // If any files need to be zero-filled, do so now
           for(;i_file_write<n_files_out;i_file_write++){
              if(n_files_out==0 && i_file_write>0)
                 SID_trap_error("Trying to create a second file in a non-multifile dataset.",ERROR_LOGIC);
              if(n_files_out>1){
                 if(i_file_write==0)
                    mkdir(filename_items_out,02755);
                 sprintf(filename_out,"%s/%s.%d",filename_items_out,filename_items_out_base,i_file_write);
              }
              else
                 sprintf(filename_out,"%s",filename_items_out);
              if(i_file_write==(n_files_out-1))
                 n_halos_write=n_halos_all-i_halo_write;
              else
                 n_halos_write=(int)((float)(i_file_write+1)*(float)n_halos_all/(float)n_files_out)-j_halo_write;
              int n_halos_left=n_halos_all-j_halo_write;
              if(n_halos_write>n_halos_left)
                 n_halos_write=n_halos_left;
              if(fp_write!=NULL)
                 fclose(fp_write);
              if((fp_write=fopen(filename_out,"w"))!=NULL){
                 fwrite(&i_file_write, sizeof(int),1,fp_write);
                 fwrite(&n_files_out,  sizeof(int),1,fp_write);
                 fwrite(&n_halos_write,sizeof(int),1,fp_write);
                 fwrite(&n_halos_all,  sizeof(int),1,fp_write);
              }
              else
                 SID_trap_error("Could not create file {%s}.",ERROR_IO_OPEN,filename_out);
           }
           if(fp_write!=NULL)
              fclose(fp_write);
           SID_log("Done.",SID_LOG_CLOSE);
        }
        SID_log("Done.",SID_LOG_CLOSE);
     } // i_run
     SID_log("Done.",SID_LOG_CLOSE);
  } // i_snap
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}

