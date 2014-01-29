#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void read_forests(char  *filename_root_in,
                  int    n_trees_group,
                  int    n_trees_subgroup,
                  int   *n_forests_group,
                  int   *n_forests_subgroup,
                  int  **i_forest_group,
                  int  **i_forest_subgroup,
                  int  **n_halos_forest_group,
                  int  **n_halos_forest_subgroup,
                  int   *n_trees_forest_groups_max,
                  int   *n_trees_forest_subgroups_max){

  int    *n_halos_tree_group;
  int    *n_halos_tree_subgroup;
  char    filename_in[MAX_FILENAME_LENGTH];
  FILE   *fp_in=NULL;
  char   *line=NULL;
  size_t  line_length=0;

  SID_log("Reading tree->forest mappings...",SID_LOG_OPEN);

  // Read mapping for groups
  int i_tree;
  int n_trees_group_in;
  sprintf(filename_in,"%s/tree2forest_mapping_groups.txt",filename_root_in);
  fp_in=fopen(filename_in,"r");
  n_trees_group_in=count_lines_data(fp_in);
  if(n_trees_group!=n_trees_group_in)
     SID_trap_error("The group tree count does not match the input (%d!=%d)",ERROR_LOGIC,n_trees_group,n_trees_group_in);
  n_halos_tree_group=(int *)SID_malloc(sizeof(int)*n_trees_group);
  (*i_forest_group) =(int *)SID_malloc(sizeof(int)*n_trees_group);
  for(i_tree=0;i_tree<n_trees_group;i_tree++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_int(line,2,&((*i_forest_group)[i_tree]));
     grab_int(line,3,&(n_halos_tree_group[i_tree]));
  }
  fclose(fp_in);

  // Count the number of group forests
  int     i_forest;
  int     i_tree_start;
  int     n_trees_forest;
  int     i_forest_group_largest;
  // Sort the forest IDs
  size_t *i_forest_group_index=NULL;
  (*n_trees_forest_groups_max)=0;
  merge_sort((*i_forest_group),(size_t)n_trees_group,&i_forest_group_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  // Skip trees with forest IDs < 0
  i_tree=0;
  while((*i_forest_group)[i_forest_group_index[i_tree]]<0 && i_tree<(n_trees_group-1))
     i_tree++;
  if((*i_forest_group)[i_forest_group_index[i_tree]]<0) i_tree++;
  for(i_forest=0;i_tree<n_trees_group;i_forest++){
     // Count how many trees are in this forest
     i_tree_start=i_tree;
     while((*i_forest_group)[i_forest_group_index[i_tree]]==i_forest && i_tree<(n_trees_group-1))
        i_tree++;
     if((*i_forest_group)[i_forest_group_index[i_tree]]==i_forest) i_tree++;
     n_trees_forest=i_tree-i_tree_start;
     if(n_trees_forest>(*n_trees_forest_groups_max)){
        (*n_trees_forest_groups_max)=n_trees_forest;
        i_forest_group_largest=i_forest;
     }
     if(n_trees_forest<1)
        SID_trap_error("There is a gap in the group tree->forest mapping (%d is missing).",ERROR_LOGIC,i_forest);
  }
  (*n_forests_group)=i_forest;
  SID_free(SID_FARG i_forest_group_index);
  SID_log("No. of    group forests                 = %d",SID_LOG_COMMENT,(*n_forests_group));

  // Read mapping for subgroups
  int n_trees_subgroup_in;
  sprintf(filename_in,"%s/tree2forest_mapping_subgroups.txt",filename_root_in);
  fp_in=fopen(filename_in,"r");
  n_trees_subgroup_in=count_lines_data(fp_in);
  if(n_trees_subgroup!=n_trees_subgroup_in)
     SID_trap_error("The subgroup tree count does not match the input (%d!=%d)",ERROR_LOGIC,n_trees_subgroup,n_trees_subgroup_in);
  n_halos_tree_subgroup=(int *)SID_malloc(sizeof(int)*n_trees_subgroup);
  (*i_forest_subgroup) =(int *)SID_malloc(sizeof(int)*n_trees_subgroup);
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_int(line,2,&((*i_forest_subgroup)[i_tree]));
     grab_int(line,3,&(n_halos_tree_subgroup[i_tree]));
  }
  fclose(fp_in);

  // Count the number of subgroup forests
  int     i_forest_subgroup_largest;
  // Sort the forest IDs
  size_t *i_forest_subgroup_index=NULL;
  (*n_trees_forest_subgroups_max)=0;
  merge_sort((*i_forest_subgroup),(size_t)n_trees_subgroup,&i_forest_subgroup_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  // Skip trees with forest IDs < 0
  i_tree=0;
  while((*i_forest_subgroup)[i_forest_subgroup_index[i_tree]]<0 && i_tree<(n_trees_subgroup-1))
     i_tree++;
  if((*i_forest_subgroup)[i_forest_subgroup_index[i_tree]]<0) i_tree++;
  for(i_forest=0;i_tree<n_trees_subgroup;i_forest++){
     // Count how many trees are in this forest
     i_tree_start=i_tree;
     while((*i_forest_subgroup)[i_forest_subgroup_index[i_tree]]==i_forest && i_tree<(n_trees_subgroup-1))
        i_tree++;
     if((*i_forest_subgroup)[i_forest_subgroup_index[i_tree]]==i_forest) i_tree++;
     n_trees_forest=i_tree-i_tree_start;
     if(n_trees_forest>(*n_trees_forest_subgroups_max)){
        (*n_trees_forest_subgroups_max)=n_trees_forest;
        i_forest_subgroup_largest=i_forest;
     }
     if(n_trees_forest<1)
        SID_trap_error("There is a gap in the subgroup tree->forest mapping (%d is missing).",ERROR_LOGIC,i_forest);
  }
  (*n_forests_subgroup)=i_forest;
  SID_free(SID_FARG i_forest_subgroup_index);
  SID_log("No. of subgroup forests                 = %d",SID_LOG_COMMENT,(*n_forests_subgroup));

  // Turn the tree halo counts into forest halo counts
  (*n_halos_forest_group)   =(int *)SID_calloc(sizeof(int)*(*n_forests_group));
  (*n_halos_forest_subgroup)=(int *)SID_calloc(sizeof(int)*(*n_forests_subgroup));
  for(i_tree=0;i_tree<n_trees_group;i_tree++)
     (*n_halos_forest_group)[(*i_forest_group)[i_tree]]+=n_halos_tree_group[i_tree];
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++)
     (*n_halos_forest_subgroup)[(*i_forest_subgroup)[i_tree]]+=n_halos_tree_subgroup[i_tree];
  SID_log("No. of trees in largest group    forest = %4d (%d halos)",SID_LOG_COMMENT,
                                                                     (*n_trees_forest_groups_max),
                                                                     (*n_halos_forest_group)[i_forest_group_largest]);
  SID_log("No. of trees in largest subgroup forest = %4d (%d halos)",SID_LOG_COMMENT,
                                                                     (*n_trees_forest_subgroups_max),
                                                                     (*n_halos_forest_subgroup)[i_forest_subgroup_largest]);

  // Clean-up
  SID_free(SID_FARG n_halos_tree_group);
  SID_free(SID_FARG n_halos_tree_subgroup);

  SID_log("Done.",SID_LOG_CLOSE);
}

