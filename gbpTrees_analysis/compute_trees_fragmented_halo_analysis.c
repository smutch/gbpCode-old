#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void compute_trees_fragmented_halo_analysis(tree_info *trees,char *filename_out_root){

  // Compute merger rates ...
  SID_log("Performing fragmented halo analysis...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Count the number of fragmented halos and the number of branches with fragmented halos
  SID_log("Counting fragmented halos and branches...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_fragmented                     =0;
  int n_fragmented_central             =0;
  int n_fragmented_substructure        =0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_fragmented(current_halo)){
           if(check_treenode_if_central(current_halo))
              n_fragmented_central++;
           else
              n_fragmented_substructure++;
           n_fragmented++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  char filename_root_branches[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,   "%s_fragmented_halos",   filename_out_root);
  sprintf(filename_root_branches,"%s_fragmented_branches",filename_out_root);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  treenode_list_info *list_halos_central;
  treenode_list_info *list_halos_substructure;
  int                *branch_number;
  double             *delta_fragmented;
  double             *Mvir_parent;
  init_treenode_list("fragmented_all",         n_fragmented,             &list_halos);
  init_treenode_list("fragmented_centrals",    n_fragmented_central,     &list_halos_central);
  init_treenode_list("fragmented_substructure",n_fragmented_substructure,&list_halos_substructure);
  init_treenode_info_data(list_halos,SID_FARG branch_number,   SID_INT,   "branch_number");
  init_treenode_info_data(list_halos,SID_FARG delta_fragmented,SID_DOUBLE,"delta_fragmented [Gyrs]");
  init_treenode_info_data(list_halos,SID_FARG Mvir_parent,     SID_DOUBLE,"M_parent [h^{-1] M_sol]");

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_fragmented(current_halo)){
           tree_markers_info markers;
           find_treenode_markers(trees,current_halo,&markers);
           if(check_treenode_if_central(current_halo))
              add_to_treenode_list(list_halos_central,current_halo);
           else
              add_to_treenode_list(list_halos_substructure,current_halo);
           add_to_treenode_list(list_halos,current_halo);
           delta_fragmented[i_list]=fetch_treenode_delta_t(trees,current_halo,markers.branch_leaf)/S_PER_GYR;
           Mvir_parent[i_list]     =fetch_treenode_Mvir   (trees,current_halo->parent);
           i_list++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write the files
  write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_markers   (trees,filename_out_root,list_halos);
  write_treenode_list_data      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_central);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_substructure);

  // Clean-up
  free_treenode_list(&list_halos);
  free_treenode_list(&list_halos_central);
  free_treenode_list(&list_halos_substructure);

  SID_log("Done.",SID_LOG_CLOSE);
}

