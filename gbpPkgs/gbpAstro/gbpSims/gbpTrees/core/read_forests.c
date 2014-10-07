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

void read_forests(const char  *filename_root_in,
                  int         *n_forests_group,
                  int         *n_forests_subgroup,
                  int         *n_forests_group_local,
                  int         *n_forests_subgroup_local,
                  int        **tree2forest_mapping_group,
                  int        **tree2forest_mapping_subgroup,
                  int         *n_trees_forest_groups_max,
                  int         *n_trees_forest_subgroups_max,
                  int         *forest_lo_group_local,
                  int         *forest_hi_group_local,
                  int         *forest_lo_subgroup_local,
                  int         *forest_hi_subgroup_local,
                  int         *n_groups_local,
                  int         *n_subgroups_local,
                  int         *n_groups_max_snap_local,
                  int         *n_subgroups_max_snap_local){

  int    *n_halos_tree_group;
  int    *n_halos_tree_subgroup;
  char    filename_in[MAX_FILENAME_LENGTH];
  FILE   *fp_in      =NULL;
  char   *line       =NULL;
  size_t  line_length=0;

  SID_log("Reading tree->forest mappings...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Read mapping for groups
  int i_tree;
  int n_trees_group;
  int n_groups_total=0;
  sprintf(filename_in,"%s/tree2forest_mapping_groups.txt",filename_root_in);
  fp_in                       =fopen(filename_in,"r");
  n_trees_group               =count_lines_data(fp_in);
  n_halos_tree_group          =(int *)SID_malloc(sizeof(int)*n_trees_group);
  (*tree2forest_mapping_group)=(int *)SID_malloc(sizeof(int)*n_trees_group);
  for(i_tree=0;i_tree<n_trees_group;i_tree++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_int(line,2,&((*tree2forest_mapping_group)[i_tree]));
     grab_int(line,3,&(n_halos_tree_group[i_tree]));
     n_groups_total+=n_halos_tree_group[i_tree];
  }
  fclose(fp_in);

  // Count the number of group forests
  int     i_forest;
  int     i_tree_start;
  int     n_trees_forest;
  int     tree2forest_mapping_group_largest;
  // Sort the forest IDs
  size_t *tree2forest_mapping_group_index=NULL;
  (*n_trees_forest_groups_max)=0;
  merge_sort((*tree2forest_mapping_group),(size_t)n_trees_group,&tree2forest_mapping_group_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  // Skip trees with forest IDs < 0
  i_tree=0;
  while((*tree2forest_mapping_group)[tree2forest_mapping_group_index[i_tree]]<0 && i_tree<(n_trees_group-1))
     i_tree++;
  if((*tree2forest_mapping_group)[tree2forest_mapping_group_index[i_tree]]<0) i_tree++;
  for(i_forest=0;i_tree<n_trees_group;i_forest++){
     // Count how many trees are in this forest
     i_tree_start=i_tree;
     while((*tree2forest_mapping_group)[tree2forest_mapping_group_index[i_tree]]==i_forest && i_tree<(n_trees_group-1))
        i_tree++;
     if((*tree2forest_mapping_group)[tree2forest_mapping_group_index[i_tree]]==i_forest) i_tree++;
     n_trees_forest=i_tree-i_tree_start;
     if(n_trees_forest>(*n_trees_forest_groups_max)){
        (*n_trees_forest_groups_max)=n_trees_forest;
        tree2forest_mapping_group_largest=i_forest;
     }
     if(n_trees_forest<1)
        SID_trap_error("There is a gap in the group tree->forest mapping (%d is missing).",ERROR_LOGIC,i_forest);
  }
  (*n_forests_group)=i_forest;
  SID_free(SID_FARG tree2forest_mapping_group_index);
  SID_log("No. of    group forests                 = %d (%d halos)",SID_LOG_COMMENT,(*n_forests_group),n_groups_total);

  // Read mapping for subgroups
  int n_trees_subgroup;
  int n_subgroups_total=0;
  sprintf(filename_in,"%s/tree2forest_mapping_subgroups.txt",filename_root_in);
  fp_in                          =fopen(filename_in,"r");
  n_trees_subgroup               =count_lines_data(fp_in);
  n_halos_tree_subgroup          =(int *)SID_malloc(sizeof(int)*n_trees_subgroup);
  (*tree2forest_mapping_subgroup)=(int *)SID_malloc(sizeof(int)*n_trees_subgroup);
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
     grab_next_line_data(fp_in,&line,&line_length);
     grab_int(line,2,&((*tree2forest_mapping_subgroup)[i_tree]));
     grab_int(line,3,&(n_halos_tree_subgroup[i_tree]));
     n_subgroups_total+=n_halos_tree_subgroup[i_tree];
  }
  fclose(fp_in);

  // Count the number of subgroup forests
  int     tree2forest_mapping_subgroup_largest;
  // Sort the forest IDs
  size_t *tree2forest_mapping_subgroup_index=NULL;
  (*n_trees_forest_subgroups_max)=0;
  merge_sort((*tree2forest_mapping_subgroup),(size_t)n_trees_subgroup,&tree2forest_mapping_subgroup_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  // Skip trees with forest IDs < 0
  i_tree=0;
  while((*tree2forest_mapping_subgroup)[tree2forest_mapping_subgroup_index[i_tree]]<0 && i_tree<(n_trees_subgroup-1))
     i_tree++;
  if((*tree2forest_mapping_subgroup)[tree2forest_mapping_subgroup_index[i_tree]]<0) i_tree++;
  for(i_forest=0;i_tree<n_trees_subgroup;i_forest++){
     // Count how many trees are in this forest
     i_tree_start=i_tree;
     while((*tree2forest_mapping_subgroup)[tree2forest_mapping_subgroup_index[i_tree]]==i_forest && i_tree<(n_trees_subgroup-1))
        i_tree++;
     if((*tree2forest_mapping_subgroup)[tree2forest_mapping_subgroup_index[i_tree]]==i_forest) i_tree++;
     n_trees_forest=i_tree-i_tree_start;
     if(n_trees_forest>(*n_trees_forest_subgroups_max)){
        (*n_trees_forest_subgroups_max)=n_trees_forest;
        tree2forest_mapping_subgroup_largest=i_forest;
     }
     if(n_trees_forest<1)
        SID_trap_error("There is a gap in the subgroup tree->forest mapping (%d is missing).",ERROR_LOGIC,i_forest);
  }
  (*n_forests_subgroup)=i_forest;
  SID_free(SID_FARG tree2forest_mapping_subgroup_index);
  SID_log("No. of subgroup forests                 = %d (%d halos)",SID_LOG_COMMENT,(*n_forests_subgroup),n_subgroups_total);

  // Turn the tree halo counts into forest halo counts
  int *n_halos_forest_group   =(int *)SID_calloc(sizeof(int)*(*n_forests_group));
  int *n_halos_forest_subgroup=(int *)SID_calloc(sizeof(int)*(*n_forests_subgroup));
  for(i_tree=0;i_tree<n_trees_group;i_tree++){
     if((*tree2forest_mapping_group)[i_tree]>=0)
        n_halos_forest_group[(*tree2forest_mapping_group)[i_tree]]+=n_halos_tree_group[i_tree];
  }
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
     if((*tree2forest_mapping_subgroup)[i_tree]>=0)
        n_halos_forest_subgroup[(*tree2forest_mapping_subgroup)[i_tree]]+=n_halos_tree_subgroup[i_tree];
  }
  SID_log("No. of trees in largest group    forest = %d (%d halos)",SID_LOG_COMMENT,
                                                                    (*n_trees_forest_groups_max),
                                                                    n_halos_forest_group[tree2forest_mapping_group_largest]);
  SID_log("No. of trees in largest subgroup forest = %d (%d halos)",SID_LOG_COMMENT,
                                                                    (*n_trees_forest_subgroups_max),
                                                                    n_halos_forest_subgroup[tree2forest_mapping_subgroup_largest]);

  // Clean-up
  SID_free(SID_FARG n_halos_tree_group);
  SID_free(SID_FARG n_halos_tree_subgroup);

  // Perform domain decomosition - groups
  int  tree_count_group_local;
  int *tree_count_split=NULL;
  int *forest_lo_split =NULL;
  int *forest_hi_split =NULL;
  split_forests_n_ways(n_halos_forest_group,
                       (*n_forests_group),
                       SID.n_proc,
                       &tree_count_split,
                       &forest_lo_split,
                       &forest_hi_split);
  tree_count_group_local  =tree_count_split[SID.My_rank];
  (*forest_lo_group_local)=forest_lo_split[SID.My_rank];
  (*forest_hi_group_local)=forest_hi_split[SID.My_rank];
  SID_free(SID_FARG tree_count_split);
  SID_free(SID_FARG forest_lo_split);
  SID_free(SID_FARG forest_hi_split);

  // Perform domain decomosition - subgroups
  int tree_count_subgroup_local;
  split_forests_n_ways(n_halos_forest_subgroup,
                       (*n_forests_subgroup),
                       SID.n_proc,
                       &tree_count_split,
                       &forest_lo_split,
                       &forest_hi_split);
  tree_count_subgroup_local  =tree_count_split[SID.My_rank];
  (*forest_lo_subgroup_local)=forest_lo_split[SID.My_rank];
  (*forest_hi_subgroup_local)=forest_hi_split[SID.My_rank];
  SID_free(SID_FARG tree_count_split);
  SID_free(SID_FARG forest_lo_split);
  SID_free(SID_FARG forest_hi_split);

  // Determine local forest information - groups
  int n_forests_group_in;
  sprintf(filename_in,"%s/forest_info_groups.txt",filename_root_in);
  fp_in=fopen(filename_in,"r");
  n_forests_group_in=count_lines_data(fp_in);
  if((*n_forests_group)!=n_forests_group_in)
     SID_trap_error("The group forest count does not match the input (%d!=%d)",ERROR_LOGIC,n_forests_group,n_forests_group_in);
  (*n_groups_local)         =0;
  (*n_groups_max_snap_local)=0;
  for(i_forest=0;i_forest<(*n_forests_group);i_forest++){
     grab_next_line_data(fp_in,&line,&line_length);
     if(i_forest>=(*forest_lo_group_local) && i_forest<=(*forest_hi_group_local)){
        int n_halos_forest;
        int n_halos_forest_max_contemp;
        grab_int(line,3,&n_halos_forest);
        grab_int(line,4,&n_halos_forest_max_contemp);
        (*n_groups_local)         +=n_halos_forest;
        (*n_groups_max_snap_local)+=n_halos_forest_max_contemp;
     }
  }
  fclose(fp_in);

  // Determine local forest information - subgroups
  int n_forests_subgroup_in;
  sprintf(filename_in,"%s/forest_info_subgroups.txt",filename_root_in);
  fp_in=fopen(filename_in,"r");
  n_forests_subgroup_in=count_lines_data(fp_in);
  if((*n_forests_subgroup)!=n_forests_subgroup_in)
     SID_trap_error("The subgroup forest count does not match the input (%d!=%d)",ERROR_LOGIC,n_forests_subgroup,n_forests_subgroup_in);
  (*n_subgroups_local)         =0;
  (*n_subgroups_max_snap_local)=0;
  for(i_forest=0;i_forest<(*n_forests_subgroup);i_forest++){
     grab_next_line_data(fp_in,&line,&line_length);
     if(i_forest>=(*forest_lo_subgroup_local) && i_forest<=(*forest_hi_subgroup_local)){
        int n_halos_forest;
        int n_halos_forest_max_contemp;
        grab_int(line,3,&n_halos_forest);
        grab_int(line,4,&n_halos_forest_max_contemp);
        (*n_subgroups_local)         +=n_halos_forest;
        (*n_subgroups_max_snap_local)+=n_halos_forest_max_contemp;
     }
  }
  fclose(fp_in);

  // Report domain decomposition results
  if(SID.n_proc>1){
     int i_rank;
     for(i_rank=0;i_rank<SID.n_proc;i_rank++){
        int forest_lo_rank;
        int forest_hi_rank;
        int tree_count_rank;
        int n_halos_rank;
        int n_alloc_rank;
        forest_lo_rank =(*forest_lo_group_local);
        forest_hi_rank =(*forest_hi_group_local);
        tree_count_rank=tree_count_group_local;
        n_halos_rank   =(*n_groups_local);
        n_alloc_rank   =(*n_groups_max_snap_local);
        SID_Bcast(&forest_lo_rank, sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&forest_hi_rank, sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&tree_count_rank,sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&n_halos_rank,   sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&n_alloc_rank,   sizeof(int),i_rank,SID.COMM_WORLD);
        SID_log("Rank #%4d will process %5d    group forests (%5d trees, %8d halos in total; alloc for %8d halos)",SID_LOG_COMMENT,
                i_rank,forest_hi_rank-forest_lo_rank+1,tree_count_rank,n_halos_rank,n_alloc_rank);
     }
     for(i_rank=0;i_rank<SID.n_proc;i_rank++){
        int forest_lo_rank;
        int forest_hi_rank;
        int tree_count_rank;
        int n_halos_rank;
        int n_alloc_rank;
        forest_lo_rank =(*forest_lo_subgroup_local);
        forest_hi_rank =(*forest_hi_subgroup_local);
        tree_count_rank=tree_count_subgroup_local;
        n_halos_rank   =(*n_subgroups_local);
        n_alloc_rank   =(*n_subgroups_max_snap_local);
        SID_Bcast(&forest_lo_rank, sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&forest_hi_rank, sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&tree_count_rank,sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&n_halos_rank,   sizeof(int),i_rank,SID.COMM_WORLD);
        SID_Bcast(&n_alloc_rank,   sizeof(int),i_rank,SID.COMM_WORLD);
        SID_log("Rank #%4d will process %5d subgroup forests (%5d trees, %8d halos in total; alloc for %8d halos)",SID_LOG_COMMENT,
                i_rank,forest_hi_rank-forest_lo_rank+1,tree_count_rank,n_halos_rank,n_alloc_rank);
     }
  }
  (*n_forests_group_local)   =(*forest_hi_group_local)   -(*forest_lo_group_local)   +1;
  (*n_forests_subgroup_local)=(*forest_hi_subgroup_local)-(*forest_lo_subgroup_local)+1;

  // Clean-up
  SID_free(SID_FARG n_halos_forest_group);
  SID_free(SID_FARG n_halos_forest_subgroup);
  SID_free(SID_FARG line);

  SID_log("Done.",SID_LOG_CLOSE);
}

