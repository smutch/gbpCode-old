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

void compute_trees_strayed_halo_analysis(tree_info *trees,char *filename_out_root_in,int i_type){

  char group_prefix[8];
  if(i_type==0)
     sprintf(group_prefix,"");
  else
     sprintf(group_prefix,"sub");

  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s_%sgroup",filename_out_root_in,group_prefix);

  // Compute merger rates ...
  SID_log("Performing strayed %sgroup analysis...",SID_LOG_OPEN|SID_LOG_TIMER,group_prefix);

  // Count the number of strayed halos 
  SID_log("Counting strayed halos...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_strayed                     =0;
  int n_strayed_central             =0;
  int n_strayed_substructure        =0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;
     if(i_type==0)
        current_halo=trees->first_neighbour_groups[i_snap];
     else
        current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_STRAYED)){
           if(check_treenode_if_central(current_halo))
              n_strayed_central++;
           else if(check_treenode_if_satellite(current_halo))
              n_strayed_substructure++;
           n_strayed++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,"%s_strayed_%sgroups",filename_out_root,group_prefix);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  treenode_list_info *list_halos_central;
  treenode_list_info *list_halos_substructure;
  int                *branch_number;
  double             *delta_strayed;
  double             *Mvir_parent;
  init_treenode_list("strayed_all",n_strayed,&list_halos);
  init_treenode_info_data(list_halos,SID_FARG Mvir_parent,SID_DOUBLE,"M_parent [h^{-1] M_sol]");
  if(i_type==1){
     init_treenode_list("strayed_centrals",    n_strayed_central,     &list_halos_central);
     init_treenode_list("strayed_substructure",n_strayed_substructure,&list_halos_substructure);
  }

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;
     if(i_type==0)
        current_halo=trees->first_neighbour_groups[i_snap];
     else
        current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_STRAYED)){
           if(check_treenode_if_central(current_halo))
              add_to_treenode_list(list_halos_central,current_halo);
           else if(check_treenode_if_satellite(current_halo))
              add_to_treenode_list(list_halos_substructure,current_halo);
           add_to_treenode_list(list_halos,current_halo);
           Mvir_parent[i_list]=fetch_treenode_Mvir(trees,current_halo->parent);
           i_list++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  finalize_treenode_list(trees,list_halos);
  if(i_type==1){
     finalize_treenode_list(trees,list_halos_central);
     finalize_treenode_list(trees,list_halos_substructure);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write the files
  //write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_data      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos);
  if(i_type==1){
     write_treenode_list_hist(trees,filename_out_root,list_halos_central);
     write_treenode_list_hist(trees,filename_out_root,list_halos_substructure);
  }

  // Clean-up
  free_treenode_list(&list_halos);
  if(i_type==1){
     free_treenode_list(&list_halos_central);
     free_treenode_list(&list_halos_substructure);
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

