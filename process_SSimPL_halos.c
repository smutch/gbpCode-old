#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpRender.h>

void process_SSimPL_halos(char *filename_SSimPL_root,
                          char *filename_halos_version,
                          char *filename_trees_version,
                          int   i_read,
                          int   mode,
                          int   select_function(int                i_group,
                                                int                j_subgroup,
                                                int                i_subgroup,
                                                int                flag_long_ids,
                                                process_halo_info *group_i,
                                                process_halo_info *subgroup_i,
                                                void              *params),
                          int   action_function(int                i_group,
                                                int                j_subgroup,
                                                int                i_subgroup,
                                                int                flag_long_ids,
                                                process_halo_info *group_i,
                                                process_halo_info *subgroup_i,
                                                void              *params),
                          void *params){

  // Set the halo and tree filename roots
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_root,filename_trees_version);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_root,filename_halos_version);

  // Initialize filename paths
  char filename_input_file_root[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_groups[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_subgroups[MAX_FILENAME_LENGTH];
  char filename_input_dir_horizontal_trees[MAX_FILENAME_LENGTH];
  strcpy(filename_input_file_root,filename_trees_root);
  strip_path(filename_input_file_root);
  sprintf(filename_input_dir_horizontal,          "%s/horizontal",filename_trees_root);
  sprintf(filename_input_dir_horizontal_groups,   "%s/groups",    filename_input_dir_horizontal);
  sprintf(filename_input_dir_horizontal_subgroups,"%s/subgroups", filename_input_dir_horizontal);
  sprintf(filename_input_dir_horizontal_trees,    "%s/trees",     filename_input_dir_horizontal);

  // Open horizontal tree file
  char filename_in[MAX_FILENAME_LENGTH];
  SID_fp fp_trees_in;
  sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_input_dir_horizontal_trees,i_read);
  SID_fopen(filename_in,"r",&fp_trees_in);

  // Open the halo files 
  char   filename_input_halos_groups[MAX_FILENAME_LENGTH];
  char   filename_input_halos_subgroups[MAX_FILENAME_LENGTH];
  char   filename_ids[MAX_FILENAME_LENGTH];
  int    n_groups_cat;
  int    n_subgroups_cat;
  int    offset_size;
  int    id_byte_size;
  SID_fp fp_groups_length;
  SID_fp fp_groups_offset;
  SID_fp fp_subgroups_length;
  SID_fp fp_subgroups_offset;
  sprintf(filename_input_halos_groups,   "%s_%03d.catalog_groups",   filename_halos_root,i_read);
  sprintf(filename_input_halos_subgroups,"%s_%03d.catalog_subgroups",filename_halos_root,i_read);
  sprintf(filename_ids,                  "%s_%03d.catalog_particles",filename_halos_root,i_read);
  SID_fopen(filename_input_halos_groups,   "r",&fp_groups_length);
  SID_fopen(filename_input_halos_groups,   "r",&fp_groups_offset);
  SID_fopen(filename_input_halos_subgroups,"r",&fp_subgroups_length);
  SID_fopen(filename_input_halos_subgroups,"r",&fp_subgroups_offset);
  SID_fread_all(&n_groups_cat,   sizeof(int),1,&fp_groups_length);
  SID_fread_all(&offset_size,    sizeof(int),1,&fp_groups_length);
  SID_fread_all(&n_subgroups_cat,sizeof(int),1,&fp_subgroups_length);
  SID_fread_all(&offset_size,    sizeof(int),1,&fp_subgroups_length);
  SID_fskip(sizeof(int),2+n_groups_cat,   &fp_groups_offset);
  SID_fskip(sizeof(int),2+n_subgroups_cat,&fp_subgroups_offset);

  // Read tree file header
  int n_step_in;
  int n_search_in;
  int n_groups;
  int n_subgroups;
  int n_groups_max_in;
  int n_subgroups_max_in;
  int n_trees_subgroup;
  int n_trees_group;
  int n_progenitors_max;
  SID_fread_all(&n_step_in,         sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_search_in,       sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_groups,          sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_trees_in);
  SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_trees_in);

  // Initialize read buffers
  SID_fp_buffer *fp_subgroups_length_buffer=NULL;
  SID_fp_buffer *fp_subgroups_offset_buffer=NULL;
  SID_fp_buffer *fp_groups_length_buffer   =NULL;
  SID_fp_buffer *fp_groups_offset_buffer   =NULL;
  SID_fp_buffer *fp_trees_in_buffer        =NULL;
  size_t n_bytes_trees                     =sizeof(int)*((7*(size_t)n_groups)+(6*(size_t)n_subgroups));
  init_SID_fp_buffer(&fp_subgroups_length,(size_t)n_subgroups*sizeof(int),SIZE_OF_MEGABYTE,&fp_subgroups_length_buffer);
  init_SID_fp_buffer(&fp_subgroups_offset,(size_t)n_subgroups*offset_size,SIZE_OF_MEGABYTE,&fp_subgroups_offset_buffer);
  init_SID_fp_buffer(&fp_groups_length,   (size_t)n_groups   *sizeof(int),SIZE_OF_MEGABYTE,&fp_groups_length_buffer);
  init_SID_fp_buffer(&fp_groups_offset,   (size_t)n_groups   *offset_size,SIZE_OF_MEGABYTE,&fp_groups_offset_buffer);
  init_SID_fp_buffer(&fp_trees_in,        n_bytes_trees,                  SIZE_OF_MEGABYTE,&fp_trees_in_buffer);

  // Open IDs file and initialize the IDs array
  int     flag_long_ids=TRUE;
  size_t  n_ids        =0;
  SID_fp  fp_ids;
  if(mode){
     SID_fopen(filename_ids,"r",&fp_ids);
     SID_fread_all(&id_byte_size,sizeof(int),1,&fp_ids);
     if(id_byte_size==sizeof(int)){
        int n_ids_i;
        SID_fread_all(&n_ids_i,sizeof(int),1,&fp_ids);
        n_ids=(size_t)n_ids_i;
        flag_long_ids=FALSE;
     }
     else
        SID_fread_all(&n_ids,sizeof(size_t),1,&fp_ids);
  }
  // Find size of largest group
  int largest_group=0;
  for(int i_group=0;i_group<n_groups;i_group++){
    int n_particles_group;SID_fread_all_buffer(&n_particles_group,sizeof(int),1,fp_groups_length_buffer);
    largest_group=MAX(largest_group,n_particles_group);
  }

  // Allocate IDs
  void *ids=NULL;
  if(mode)
     ids=SID_malloc(id_byte_size*largest_group);
  int    *ids_i=(int    *)ids;
  size_t *ids_l=(size_t *)ids;
  // Reset the group length pointer
  SID_fseek(&fp_groups_length,sizeof(int),2,SID_SEEK_SET);
  reset_SID_fp_buffer(&fp_groups_length_buffer);

  // Initialize the datastructures which will hold the group and subgroup information
  process_halo_info group_i;
  process_halo_info subgroup_i;
  group_i.ids           =ids;
  subgroup_i.n_subgroups=0;
  
  // Read each group in turn
  //pcounter_info pcounter;
  //SID_init_pcounter(&pcounter,n_groups,10);
  size_t index_last_read        =0;
  int    n_groups_unused        =0;
  int    n_groups_added_multiply=0;
  int    i_group;
  int    i_subgroup;
  for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
    // Read horizontal trees for groups
    SID_fread_all_buffer(&(group_i.n_particles),  sizeof(int),1,fp_groups_length_buffer);
    SID_fread_all_buffer(&(group_i.id),           sizeof(int),1,fp_trees_in_buffer);
    SID_fread_all_buffer(&(group_i.tree_case),    sizeof(int),1,fp_trees_in_buffer);
    SID_fread_all_buffer(&(group_i.descendant_id),sizeof(int),1,fp_trees_in_buffer);
    SID_fread_all_buffer(&(group_i.tree_id),      sizeof(int),1,fp_trees_in_buffer);
    SID_fread_all_buffer(&(group_i.file_offset),  sizeof(int),1,fp_trees_in_buffer);
    SID_fread_all_buffer(&(group_i.file_index),   sizeof(int),1,fp_trees_in_buffer);
    SID_fread_all_buffer(&(group_i.n_subgroups),  sizeof(int),1,fp_trees_in_buffer);
    switch(flag_long_ids){
       case TRUE:
          SID_fread_all_buffer(&(group_i.offset),offset_size,1,fp_groups_offset_buffer);
          break;
       case FALSE:
          int offset_i;
          SID_fread_all_buffer(&offset_i,offset_size,1,fp_groups_offset_buffer);
          group_i.offset=(size_t)offset_i;
          break;
    }

    // Read each subgroup in turn
    int flag_group_ids_unread=TRUE;
    int j_subgroup;
    int i_subgroup_valid;
    tree_node_info *group_node;
    tree_node_info *subgroup_node;
    for(j_subgroup=0;j_subgroup<group_i.n_subgroups;i_subgroup++,j_subgroup++){
       // Read horizontal trees for subgroups
       SID_fread_all_buffer(&(subgroup_i.n_particles),  sizeof(int),1,fp_subgroups_length_buffer);
       SID_fread_all_buffer(&(subgroup_i.id),           sizeof(int),1,fp_trees_in_buffer);
       SID_fread_all_buffer(&(subgroup_i.tree_case),    sizeof(int),1,fp_trees_in_buffer);
       SID_fread_all_buffer(&(subgroup_i.descendant_id),sizeof(int),1,fp_trees_in_buffer);
       SID_fread_all_buffer(&(subgroup_i.tree_id),      sizeof(int),1,fp_trees_in_buffer);
       SID_fread_all_buffer(&(subgroup_i.file_offset),  sizeof(int),1,fp_trees_in_buffer);
       SID_fread_all_buffer(&(subgroup_i.file_index),   sizeof(int),1,fp_trees_in_buffer);
       switch(flag_long_ids){
          case TRUE:
             SID_fread_all_buffer(&(subgroup_i.offset),offset_size,1,fp_subgroups_offset_buffer);
             break;
          case FALSE:
             int offset_i;
             SID_fread_all_buffer(&offset_i,offset_size,1,fp_subgroups_offset_buffer);
             subgroup_i.offset=(size_t)offset_i;
             break;
       }
       // Perform selection and action
       if(select_function(i_group,
                          j_subgroup,
                          i_subgroup,
                          flag_long_ids,
                          &group_i,
                          &subgroup_i,
                          params)){
          // Read IDs
          if(mode){
             if(flag_group_ids_unread){
                SID_fseek(&fp_ids,(off_t)(group_i.offset-index_last_read),id_byte_size,SEEK_CUR);
                SID_fread_all(ids,id_byte_size,group_i.n_particles,&fp_ids);
                index_last_read=(size_t)group_i.offset+(size_t)group_i.n_particles;
                flag_group_ids_unread=FALSE;
             }
             switch(flag_long_ids){
                case TRUE:
                   subgroup_i.ids=&(ids_l[subgroup_i.offset-group_i.offset]);
                   break;
                case FALSE:
                   subgroup_i.ids=&(ids_i[subgroup_i.offset-group_i.offset]);
                   break;
             }
          }
          action_function(i_group,
                          j_subgroup,
                          i_subgroup,
                          flag_long_ids,
                          &group_i,
                          &subgroup_i,
                          params);
       }
    }
    //SID_check_pcounter(&pcounter,i_group);
  } // i_group

  // Free the buffers and perform sanity checks
  SID_fclose(&fp_trees_in);
  SID_fclose(&fp_groups_length);
  SID_fclose(&fp_groups_offset);
  SID_fclose(&fp_subgroups_length);
  SID_fclose(&fp_subgroups_offset);
  free_SID_fp_buffer(&fp_trees_in_buffer);
  free_SID_fp_buffer(&fp_groups_length_buffer);
  free_SID_fp_buffer(&fp_groups_offset_buffer);
  free_SID_fp_buffer(&fp_subgroups_length_buffer);
  free_SID_fp_buffer(&fp_subgroups_offset_buffer);
  if(mode){
     SID_fclose(&fp_ids);
     SID_free(SID_FARG ids);
  }
}

