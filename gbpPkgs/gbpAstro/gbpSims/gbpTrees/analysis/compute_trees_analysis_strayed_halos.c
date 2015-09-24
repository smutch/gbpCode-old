#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void compute_trees_analysis_strayed_halos(tree_info *trees,char *filename_out_root_in,int i_type,double logM_min,double dlogM,int n_logM){

  tree_node_info **neighbour_list_start=NULL;
  char group_prefix[8];
  if(i_type==0){
     sprintf(group_prefix,"");
     neighbour_list_start=trees->first_neighbour_groups;
  }
  else{
     sprintf(group_prefix,"sub");
     neighbour_list_start=trees->first_neighbour_subgroups;
  }

  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s_%sgroup",filename_out_root_in,group_prefix);

  SID_log("Performing strayed %sgroup analysis...",SID_LOG_OPEN|SID_LOG_TIMER,group_prefix);

  // Count the number of strayed halos and the number of branches with strayed halos
  SID_log("Counting strayed halos...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_strayed=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=neighbour_list_start[i_snap];
     while(current_halo!=NULL){
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_STRAYED) && current_halo->descendant==NULL)
           n_strayed++;
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  char filename_root_branches[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,   "%s_strayed_%sgroup",         filename_out_root,group_prefix);
  sprintf(filename_root_branches,"%s_strayed_%sgroup_branches",filename_out_root,group_prefix);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  init_treenode_list("strayed",n_strayed,&list_halos);
  double *delta_normed;   init_treenode_info_data(list_halos,SID_FARG delta_normed,   SID_DOUBLE,"Strayed branch length/t_dyn(z)");
  double *delta_strayed;  init_treenode_info_data(list_halos,SID_FARG delta_strayed,  SID_DOUBLE,"Strayed branch length [Gyrs]");
  int    *delta_snap;     init_treenode_info_data(list_halos,SID_FARG delta_snap,     SID_INT,   "Strayed branch length (snapshots)");

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=neighbour_list_start[i_snap];
     double          z           =(double)trees->z_list[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_STRAYED) && current_halo->descendant==NULL){
           // Add halo to list
           add_to_treenode_list(list_halos,current_halo);

           // Fetch halo markers
           tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,current_halo);

           // Compile strayed halo data
           delta_normed     [i_list]=fetch_treenode_delta_t      (trees,current_halo,markers->branch_leaf)/t_dyn_z(z,trees->cosmo);
           delta_strayed    [i_list]=fetch_treenode_delta_t      (trees,current_halo,markers->branch_leaf)/S_PER_GYR;
           delta_snap       [i_list]=fetch_treenode_snapshot(trees,current_halo)-fetch_treenode_snapshot(trees,markers->branch_leaf);

           i_list++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }

  // Finalize list
  finalize_treenode_list(trees,list_halos);
  SID_log("Done.",SID_LOG_CLOSE);

  // Write the files
  write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_data      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos,logM_min,dlogM,n_logM);

  // Clean-up
  free_treenode_list(&list_halos);

  SID_log("Done.",SID_LOG_CLOSE);
}

