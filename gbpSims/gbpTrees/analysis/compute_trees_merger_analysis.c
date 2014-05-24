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

void compute_trees_merger_analysis(tree_info *trees,char *filename_out_root){

  // Compute merger rates ...
  SID_log("Performing merger analysis...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Count the number of mergers 
  SID_log("Counting mergers...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_mergers                     =0;
  int n_mergers_central             =0;
  int n_mergers_substructure        =0;
  int n_mergers_nofrags             =0;
  int n_mergers_nofrags_central     =0;
  int n_mergers_nofrags_substructure=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new merger
        if(check_treenode_if_merger(current_halo)){
           if(!check_treenode_if_fragmented(current_halo)){
             if(check_treenode_if_central(current_halo))
                n_mergers_nofrags_central++;
             else
                n_mergers_nofrags_substructure++;
             n_mergers_nofrags++;
           }
           if(check_treenode_if_central(current_halo))
              n_mergers_central++;
           else
              n_mergers_substructure++;
           n_mergers++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  treenode_list_info *list_halos_central;
  treenode_list_info *list_halos_substructure;
  treenode_list_info *list_halos_nofrags;
  treenode_list_info *list_halos_nofrags_central;
  treenode_list_info *list_halos_nofrags_substructure;
  init_treenode_list("mergers_all",                 n_mergers,                     &list_halos);
  init_treenode_list("mergers_centrals",            n_mergers_central,             &list_halos_central);
  init_treenode_list("mergers_substructure",        n_mergers_substructure,        &list_halos_substructure);
  init_treenode_list("mergers_nofrags_all",         n_mergers_nofrags,             &list_halos_nofrags);
  init_treenode_list("mergers_nofrags_centrals",    n_mergers_nofrags_central,     &list_halos_nofrags_central);
  init_treenode_list("mergers_nofrags_substructure",n_mergers_nofrags_substructure,&list_halos_nofrags_substructure);

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list     =0;
  int n_emerged_i=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new merger
        if(check_treenode_if_merger(current_halo)){
           if(!check_treenode_if_fragmented(current_halo)){
              if(check_treenode_if_central(current_halo))
                 add_to_treenode_list(list_halos_nofrags_central,current_halo);
              else
                 add_to_treenode_list(list_halos_nofrags_substructure,current_halo);
              add_to_treenode_list(list_halos_nofrags,current_halo);
           }
           if(check_treenode_if_central(current_halo))
              add_to_treenode_list(list_halos_central,current_halo);
           else
              add_to_treenode_list(list_halos_substructure,current_halo);
           add_to_treenode_list(list_halos,current_halo);
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write the files
  write_treenode_list_markers   (trees,filename_out_root,list_halos);
  write_treenode_list_markers   (trees,filename_out_root,list_halos_nofrags);
  write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_properties(trees,filename_out_root,list_halos_nofrags);
  write_treenode_list_hist      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_central);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_substructure);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_nofrags);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_nofrags_central);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_nofrags_substructure);

  // Clean-up
  free_treenode_list(&list_halos);
  free_treenode_list(&list_halos_central);
  free_treenode_list(&list_halos_substructure);
  free_treenode_list(&list_halos_nofrags);
  free_treenode_list(&list_halos_nofrags_central);
  free_treenode_list(&list_halos_nofrags_substructure);

  SID_log("Done.",SID_LOG_CLOSE);
}

