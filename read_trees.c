#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void read_trees(char             *filename_trees_root,
                int               i_file_start,
                int               i_file_stop,
                int               mode,
                tree_node_info  **trees){
  char    filename_trees[256];
  FILE   *fp;
  int     i,j,k,l;
  int     i_tree;
  int     i_file;
  int     max_group_id;
  int     max_subgroup_id;
  int     n_files;
  int     i_line;
  char   *line=NULL;
  size_t  line_length=0;
  int     group_id_1;
  int     group_idx_1;
  int     group_tree_1;
  int     descendant_id_1;
  int     descendant_idx_1;
  int     subdescendant_id_1;
  int     subdescendant_idx_1;
  int     n_subgroups_1;
  int     i_subgroup;
  int     subgroup_id_1;
  int     subgroup_idx_1;
  int     subgroup_tree_1;
  int    *group_idx;
  int    *group_descendant;
  int    *group_tree;
  int    *subgroup_idx;
  int    *subgroup_descendant;
  int    *subgroup_tree;
  int     j_file;
  char    filename_root_out[256];
  int     n_lines;
  size_t  n_particles;
  int     n_trees;
  size_t  n_trees_t;
  int     n_subtrees;
  size_t  n_subtrees_t;
  int    *n_subtrees_tree;
  int     n_trees_in;
  int     n_subtrees_in;
  int     n_ids;
  int    *tree_id;
  int    *tree_length;
  int    *subtree_length;
  int    *tree_offsets;
  size_t *input_id;
  int     id_in;
  long    record_length;
  char    parm_txt[256];
  double  expansion_factor;
  float   overdensity_vir;
  int     n_p_min;
  int     n_p_max;
  double *M_vir;
  int    *n_p_vir;
  float  *x_bound;
  float  *y_bound;
  float  *z_bound;
  double *x_COM;
  double *y_COM;
  double *z_COM;
  double *x_CODen;
  double *y_CODen;
  double *z_CODen;
  double *vx_COM;
  double *vy_COM;
  double *vz_COM;
  double *j_x;
  double *j_y;
  double *j_z;
  double *lambda;
  double *c_15;
  int     flag_read_trees;
  int     flag_read_ids;
  int     flag_read_subtrees;
  int     flag_read_props;

  size_t  n_dims_group[2];
  size_t  n_dims_subgroup[2];

  SID_profile_start("read_trees",SID_PROFILE_NOTMPIENABLED);

  SID_log("Reading trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Create flags from mode
  /*
    flag_read_trees   =TRUE;
    flag_read_ids     =TRUE;
    flag_read_subtrees=FALSE;
    flag_read_props   =FALSE;
    if((check_mode_for_flag(mode,READ_treeS_ALL)          ||
    check_mode_for_flag(mode,READ_treeS_SUBtreeS))    &&
    !check_mode_for_flag(mode,READ_treeS_NOSUBtreeS))
    flag_read_subtrees=TRUE;
    if((check_mode_for_flag(mode,READ_treeS_ALL)         ||
    check_mode_for_flag(mode,READ_treeS_PROPERTIES)) &&
    !check_mode_for_flag(mode,READ_treeS_NOPROPERTIES))
    flag_read_props=TRUE;
    if((check_mode_for_flag(mode,READ_treeS_NOIDS))
    flag_read_ids=FALSE;
  */

  // Count number of group and subgroup ids
  SID_log("Counting groups and subgroups...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_file=i_file_stop,max_group_id=1,max_subgroup_id=1,n_files=0;i_file>=i_file_start;i_file--,n_files++){
    sprintf(filename_trees,"%s.trees_%d",filename_trees_root,i_file);
    fp=fopen(filename_trees,"r");
    n_lines=count_lines_data(fp);
    for(i_line=0;i_line<n_lines;i_line++){
      grab_next_line_data(fp,&line,&line_length);
      grab_int(line,3,&group_id_1);
      grab_int(line,8,&n_subgroups_1);
      max_group_id=MAX(max_group_id,group_id_1);
      for(i_subgroup=0;i_subgroup<n_subgroups_1;i_subgroup++){
        grab_int(line,11+i_subgroup*7,&subgroup_id_1);
        max_subgroup_id=MAX(max_subgroup_id,subgroup_id_1);
      }
    }
    fclose(fp);
  }
  SID_log("Done. (%d groups and %d subgroups)",SID_LOG_CLOSE,max_group_id,max_subgroup_id);

  // Read catalog file
  SID_log("Allocating arrays...",SID_LOG_OPEN|SID_LOG_TIMER);
  group_idx          =(int *)SID_malloc(sizeof(int)*n_files*max_group_id);
  group_descendant   =(int *)SID_malloc(sizeof(int)*n_files*max_group_id);
  group_tree         =(int *)SID_malloc(sizeof(int)*n_files*max_group_id);
  subgroup_idx       =(int *)SID_malloc(sizeof(int)*n_files*max_subgroup_id);
  subgroup_descendant=(int *)SID_malloc(sizeof(int)*n_files*max_subgroup_id);
  subgroup_tree      =(int *)SID_malloc(sizeof(int)*n_files*max_subgroup_id);
  SID_log("Done.",SID_LOG_CLOSE);

  // Read catalog file
  SID_log("Reading tree data...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_file=i_file_stop,j_file=0;i_file>=i_file_start;i_file--,j_file++){
    sprintf(filename_trees,"%s.trees_%d",filename_trees_root,i_file);
    fp=fopen(filename_trees,"r");
    n_lines=count_lines_data(fp);
    for(i_line=0;i_line<n_lines;i_line++){
      grab_next_line_data(fp,&line,&line_length);
      grab_int(line,1,&group_idx_1);
      grab_int(line,3,&group_id_1);
      grab_int(line,4,&descendant_id_1);
      grab_int(line,5,&group_tree_1);
      grab_int(line,6,&descendant_idx_1);
      grab_int(line,8,&n_subgroups_1);
      group_idx[INDEX_2D(max_group_id,j_file,group_id_1)]       =group_idx_1;
      group_descendant[INDEX_2D(max_group_id,j_file,group_id_1)]=descendant_id_1;
      group_tree[INDEX_2D(max_group_id,j_file,group_id_1)]      =group_tree_1;
      for(i_subgroup=0;i_subgroup<n_subgroups_1;i_subgroup++){
        grab_int(line,9 +i_subgroup*7,&subgroup_idx_1);
        grab_int(line,11+i_subgroup*7,&subgroup_id_1);
        grab_int(line,12+i_subgroup*7,&subdescendant_id_1);
        grab_int(line,13+i_subgroup*7,&subgroup_tree_1);
        grab_int(line,14+i_subgroup*7,&subdescendant_idx_1);
        subgroup_idx[INDEX_2D(max_subgroup_id,j_file,subgroup_id_1)]       =subgroup_idx_1;
        subgroup_descendant[INDEX_2D(max_subgroup_id,j_file,subgroup_id_1)]=subdescendant_id_1;
        subgroup_tree[INDEX_2D(max_subgroup_id,j_file,subgroup_id_1)]      =subgroup_tree_1;
      }
    }
    fclose(fp);
  }
  SID_log("Done.",SID_LOG_CLOSE);
  
  // Store results
  ADaPS_store((ADaPS **)trees,
              (void *)&n_files,
              "n_times",
              ADaPS_COPY,
              sizeof(int));
  ADaPS_store((ADaPS **)trees,
              (void *)&max_group_id,
              "n_groups",
              ADaPS_COPY,
              sizeof(int));
  ADaPS_store((ADaPS **)trees,
              (void *)&max_subgroup_id,
              "n_subgroups",
              ADaPS_COPY,
              sizeof(int));
  ADaPS_store((ADaPS **)trees,
              (void *)group_idx,
              "group_idx",
              ADaPS_DEFAULT);
  ADaPS_store((ADaPS **)trees,
              (void *)group_descendant,
              "group_descendant",
              ADaPS_DEFAULT);
  ADaPS_store((ADaPS **)trees,
              (void *)group_tree,
              "group_tree",
              ADaPS_DEFAULT);
  ADaPS_store((ADaPS **)trees,
              (void *)subgroup_idx,
              "subgroup_idx",
              ADaPS_DEFAULT);
  ADaPS_store((ADaPS **)trees,
              (void *)subgroup_descendant,
              "subgroup_descendant",
              ADaPS_DEFAULT);
  ADaPS_store((ADaPS **)trees,
              (void *)subgroup_tree,
              "subgroup_tree",
              ADaPS_DEFAULT);

  SID_log("Done.",SID_LOG_CLOSE);
}

