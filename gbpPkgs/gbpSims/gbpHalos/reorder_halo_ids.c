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
  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_in_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  int start_snap;
  int stop_snap;
  int step_snap;
  strcpy(filename_in_root, argv[1]);
  strcpy(filename_out_root,argv[2]);
  start_snap         =atoi(argv[3]);
  stop_snap          =atoi(argv[4]);
  step_snap          =atoi(argv[4]);

  int offset_size=sizeof(unsigned int);

  SID_log("Reordering IDs with root {%s} and writing to root {%s}...",SID_LOG_OPEN,filename_in_root,filename_out_root);
  int i_snap;
  for(i_snap=start_snap;i_snap<=stop_snap;i_snap+=step_snap){
     SID_log("Processing snapshot No. %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap);

     // Open files
     char  filename_in_base[MAX_FILENAME_LENGTH];
     char  filename_groups_in[MAX_FILENAME_LENGTH];
     char  filename_subgroups_in[MAX_FILENAME_LENGTH];
     char  filename_ids_in[MAX_FILENAME_LENGTH];
     char  filename_indices_in[MAX_FILENAME_LENGTH];
     char  filename_ids_out[MAX_FILENAME_LENGTH];
     char  group_prefix_text[8];
     FILE *fp_groups_size_in=NULL;
     FILE *fp_groups_offset_in=NULL;
     FILE *fp_groups_nsubs_in=NULL;
     FILE *fp_subgroups_size_in=NULL;
     FILE *fp_subgroups_offset_in=NULL;
     FILE *fp_ids_in=NULL;
     FILE *fp_indices_in=NULL;
     FILE *fp_ids_out=NULL;
     strcpy(filename_in_base,filename_in_root);
     strip_path(filename_in_base);
     sprintf(filename_groups_in,   "%s_%03d.catalog_groups",           filename_in_root, i_snap);
     sprintf(filename_subgroups_in,"%s_%03d.catalog_subgroups",        filename_in_root, i_snap);
     sprintf(filename_ids_in,      "%s_%03d.catalog_particles",        filename_in_root, i_snap);
     sprintf(filename_ids_out,     "%s_%03d.catalog_particles",        filename_out_root,i_snap);
     if((fp_groups_size_in=fopen(filename_groups_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_groups_in);
     if((fp_groups_offset_in=fopen(filename_groups_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_groups_in);
     if((fp_groups_nsubs_in=fopen(filename_groups_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_groups_in);
     if((fp_subgroups_size_in=fopen(filename_subgroups_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_subgroups_in);
     if((fp_subgroups_offset_in=fopen(filename_subgroups_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_subgroups_in);
     if((fp_ids_in=fopen(filename_ids_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_ids_in);
     fp_ids_out=fopen(filename_ids_out,"w");

     // Deal with the indices file(s) separately because they have a multi-file format
     int flag_indices_multifile;
     int i_file_indices=0;
     int i_file_indices_in;
     int n_files_indices;
     int n_subgroups_indices_i;
     int n_subgroups_indices_all;
     sprintf(filename_indices_in,"%s_%03d.catalog_subgroups_indices/%s_%03d.catalog_subgroups_indices.%d",filename_in_root,i_snap,filename_in_base,i_snap,i_file_indices);
     if((fp_indices_in=fopen(filename_indices_in,"r"))==NULL){
        sprintf(filename_indices_in,"%s_%03d.catalog_subgroups_indices",filename_in_root,i_snap);
        if((fp_indices_in=fopen(filename_indices_in,"r"))==NULL)
           SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_indices_in);
        else
           flag_indices_multifile=FALSE;
     }
     else
        flag_indices_multifile=TRUE;
     fread_verify(&i_file_indices_in,      sizeof(int),1,fp_indices_in);
     fread_verify(&n_files_indices,        sizeof(int),1,fp_indices_in);
     fread_verify(&n_subgroups_indices_i,  sizeof(int),1,fp_indices_in);
     fread_verify(&n_subgroups_indices_all,sizeof(int),1,fp_indices_in);

     // Read (or skip) headers
     int flag_long_ids;
     int flag_long_group_offsets;
     int flag_long_subgroup_offsets;
     int n_groups;
     int n_subgroups;
     int byte_size_ids;
     int byte_size_group_offsets;
     int byte_size_subgroup_offsets;
     int n_particles;
     fread_verify(&n_groups,                  sizeof(int),1,fp_groups_size_in);
     fread_verify(&byte_size_group_offsets,   sizeof(int),1,fp_groups_size_in);
     fread_verify(&n_subgroups,               sizeof(int),1,fp_subgroups_size_in);
     fread_verify(&byte_size_subgroup_offsets,sizeof(int),1,fp_subgroups_size_in);
     fread_verify(&byte_size_ids,             sizeof(int),1,fp_ids_in);
     if(byte_size_ids==sizeof(int)){
        SID_log("(int IDs)...",SID_LOG_CONTINUE);
        int n_particles_in;
        flag_long_ids=FALSE;
        fread_verify(&n_particles_in,sizeof(int),1,fp_ids_in);
        n_particles=(size_t)n_particles_in;
     }
     else if(byte_size_ids==sizeof(long long)){
        SID_log("(long long IDs)...",SID_LOG_CONTINUE);
        flag_long_ids=TRUE;
        fread_verify(&n_particles,sizeof(size_t),1,fp_ids_in);
     }
     else
        SID_trap_error("Invalid particle ID byte size (%d).",ERROR_LOGIC,byte_size_ids);
     if(byte_size_group_offsets==sizeof(unsigned int))
        flag_long_group_offsets=FALSE;
     else if(byte_size_group_offsets==sizeof(int64_t))
        flag_long_group_offsets=TRUE;
     else
        SID_trap_error("Invalid group offset byte size (%d).",ERROR_LOGIC,byte_size_group_offsets);
     if(byte_size_subgroup_offsets==sizeof(unsigned int))
        flag_long_subgroup_offsets=FALSE;
     else if(byte_size_subgroup_offsets==sizeof(int64_t))
        flag_long_subgroup_offsets=TRUE;
     else
        SID_trap_error("Invalid subgroup offset byte size (%d).",ERROR_LOGIC,byte_size_subgroup_offsets);
     fseeko(fp_groups_offset_in,   (off_t)(2*sizeof(int)+n_groups*sizeof(int)),                          SEEK_SET);
     fseeko(fp_groups_nsubs_in,    (off_t)(2*sizeof(int)+n_groups*(sizeof(int)+byte_size_group_offsets)),SEEK_SET);
     fseeko(fp_subgroups_offset_in,(off_t)(2*sizeof(int)+n_subgroups*sizeof(int)),                       SEEK_SET);

     // Write header for new IDs file
     fwrite(&byte_size_ids,sizeof(int),1,fp_ids_out);
     if(byte_size_ids==sizeof(int)){
        int n_particles_out;
        n_particles_out=(int)n_particles;
        fwrite(&n_particles_out,sizeof(int),1,fp_ids_out);
     }
     else if(byte_size_ids==sizeof(long long)){
        long long n_particles_out;
        n_particles_out=(long long)n_particles;
        fwrite(&n_particles_out,sizeof(long long),1,fp_ids_out);
     }

     // Loop over each group
     char      *ids_in;
     char      *ids_out;
     size_t    *indices_in;
     int        i_group;
     int        i_subgroup;
     int        k_subgroup;
     size_t     i_particle;
     int        n_particles_i;
     int        n_alloc=1024*1024;
     ids_in    =(char   *)SID_malloc(n_alloc*byte_size_ids);
     ids_out   =(char   *)SID_malloc(n_alloc*byte_size_ids);
     indices_in=(size_t *)SID_malloc(n_alloc*sizeof(size_t));
     for(i_group=0,i_particle=0,i_subgroup=0,k_subgroup=0;i_group<n_groups;i_group++,i_particle+=(size_t)n_particles_i){
	int     n_particles_i;
        int64_t offset_i;
        int     n_subs_i;

        // Read group size
        fread_verify(&n_particles_i,sizeof(int),1,fp_groups_size_in);

        // Read group offsets
        if(flag_long_group_offsets)
           fread_verify(&offset_i,byte_size_group_offsets,1,fp_groups_offset_in);
        else{
           unsigned int offset_in;
           fread_verify(&offset_in,byte_size_subgroup_offsets,1,fp_groups_offset_in);
           offset_i=(int64_t)offset_in;
        }

        // Read the number of subgroups for this group
        fread_verify(&n_subs_i,sizeof(int),1,fp_groups_nsubs_in);

        // Make sure the id and index array sizes are large enough
        while(n_particles_i>n_alloc){
           n_alloc*=2;
           ids_in     =(char   *)realloc(ids_in,    n_alloc*byte_size_ids);
           ids_out    =(char   *)realloc(ids_out,   n_alloc*byte_size_ids);
           indices_in =(size_t *)realloc(indices_in,n_alloc*sizeof(size_t));
        }

        // Read ids
        fread_verify(ids_in,byte_size_ids,n_particles_i,fp_ids_in);
        memcpy(ids_out,ids_in,byte_size_ids*n_particles_i);
       
        // Process this group's subgroups 
        int j_subgroup;
        for(j_subgroup=0;j_subgroup<n_subs_i;k_subgroup++,j_subgroup++,i_subgroup++){
           int     n_particles_sub;
           int64_t offset_sub;
           int     n_particles_sub_idx;

           // Read subgroup size
           fread_verify(&n_particles_sub,sizeof(int),1,fp_subgroups_size_in);

           // Read subgroup offset
           if(flag_long_subgroup_offsets)
              fread_verify(&offset_sub,byte_size_subgroup_offsets,1,fp_subgroups_offset_in);
           else{
              unsigned int offset_in;
              fread_verify(&offset_in,byte_size_subgroup_offsets,1,fp_subgroups_offset_in);
              offset_sub=(int64_t)offset_in;
           }
           offset_sub-=offset_i;

           // Read indices
           int n_particles_indices;
           if(i_subgroup>=n_subgroups_indices_i){
              i_file_indices++;
              if(fp_indices_in!=NULL)
                 fclose(fp_indices_in);
              sprintf(filename_indices_in,"%s_%03d.catalog_subgroups_indices/%s_%03d.catalog_subgroups_indices.%d",filename_in_root,i_snap,filename_in_base,i_snap,i_file_indices);
              if((fp_indices_in=fopen(filename_indices_in,"r"))==NULL){
                 sprintf(filename_indices_in,"%s_%03d.catalog_subgroups_indices",filename_in_root,i_snap);
                 if((fp_indices_in=fopen(filename_indices_in,"r"))==NULL)
                    SID_trap_error("Could not open {%s}",ERROR_IO_OPEN,filename_indices_in);
                 else
                    flag_indices_multifile=FALSE;
              }
              else
                 flag_indices_multifile=TRUE;
              fread_verify(&i_file_indices_in,      sizeof(int),1,fp_indices_in);
              fread_verify(&n_files_indices,        sizeof(int),1,fp_indices_in);
              fread_verify(&n_subgroups_indices_i,  sizeof(int),1,fp_indices_in);
              fread_verify(&n_subgroups_indices_all,sizeof(int),1,fp_indices_in);
              i_subgroup=0;
           }
           fread_verify(&n_particles_indices,sizeof(int),1,fp_indices_in);
           if(n_particles_indices!=n_particles_sub)
              SID_trap_error("Subgroup sizes don't match between files (ie %d!=%d) for i_group=%d/k_subgroup=%d",ERROR_LOGIC,n_particles_indices,n_particles_sub,i_group,k_subgroup);
           fread_verify(indices_in,sizeof(size_t),n_particles_sub,fp_indices_in);

           // Re-arrange subgroup indices
           size_t i_particle;
           size_t j_particle;
           if(flag_long_ids){
              long long *ids_in_temp;
              long long *ids_out_temp;
              ids_in_temp =(long long *)ids_in;
              ids_out_temp=(long long *)ids_out;
              for(j_particle=0;j_particle<n_particles_sub;j_particle++)
                 ids_out_temp[offset_sub+j_particle]=ids_in_temp[offset_sub+indices_in[j_particle]];
           }
           else{
              int *ids_in_temp;
              int *ids_out_temp;
              ids_in_temp =(int *)ids_in;
              ids_out_temp=(int *)ids_out;
              for(j_particle=0;j_particle<n_particles_sub;j_particle++)
                 ids_out_temp[offset_sub+j_particle]=ids_in_temp[offset_sub+indices_in[j_particle]];
           }
        }
        // Write to the new file
        fwrite(ids_out,byte_size_ids,n_particles_i,fp_ids_out);
     }
     SID_free(SID_FARG ids_in);
     SID_free(SID_FARG ids_out);
     SID_free(SID_FARG indices_in);

     // Close files
     fclose(fp_groups_size_in);
     fclose(fp_groups_offset_in);
     fclose(fp_groups_nsubs_in);
     fclose(fp_subgroups_size_in);
     fclose(fp_subgroups_offset_in);
     fclose(fp_ids_in);
     fclose(fp_indices_in);
     fclose(fp_ids_out);
     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}

