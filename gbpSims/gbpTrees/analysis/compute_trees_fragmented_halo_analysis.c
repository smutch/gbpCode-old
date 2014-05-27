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

void compute_trees_fragmented_halo_analysis(tree_info *trees,char *filename_out_root){

  // Compute merger rates ...
  SID_log("Performing fragmented halo analysis...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Count the number of fragmented halos 
  SID_log("Counting fragmented halos...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_fragmented                 =0;
  int n_fragmented_central         =0;
  int n_fragmented_substructure    =0;
  int n_fragmented_new             =0;
  int n_fragmented_new_central     =0;
  int n_fragmented_new_substructure=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_fragmented(current_halo)){
           int flag_new=FALSE;
           if(current_halo->n_progenitors==0)
              flag_new=TRUE;
           if(check_treenode_if_central(current_halo)){
              n_fragmented_central++;
              if(flag_new)
                 n_fragmented_new_central++;
           }
           else{
              n_fragmented_substructure++;
              if(flag_new)
                 n_fragmented_new_substructure++;
           }
           n_fragmented++;
           if(flag_new)
              n_fragmented_new++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  char filename_root_new[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,"%s_fragmented_halos",    filename_out_root);
  sprintf(filename_root_new,  "%s_fragmented_new_halos",filename_out_root);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  treenode_list_info *list_halos_central;
  treenode_list_info *list_halos_substructure;
  treenode_list_info *list_halos_new;
  treenode_list_info *list_halos_new_central;
  treenode_list_info *list_halos_new_substructure;
  int                *frag_snapshot;
  int                *frag_index;
  int                *frag_length;
  double             *delta_fragmented;
  double             *Mvir_parent;
  double             *Mvir;
  init_treenode_list("fragmented_all",             n_fragmented,                 &list_halos);
  init_treenode_list("fragmented_centrals",        n_fragmented_central,         &list_halos_central);
  init_treenode_list("fragmented_substructure",    n_fragmented_substructure,    &list_halos_substructure);
  init_treenode_list("fragmented_new_all",         n_fragmented_new,             &list_halos_new);
  init_treenode_list("fragmented_new_centrals",    n_fragmented_new_central,     &list_halos_new_central);
  init_treenode_list("fragmented_new_substructure",n_fragmented_new_substructure,&list_halos_new_substructure);
  init_treenode_info_data(list_halos_new,SID_FARG frag_length,     SID_INT,   "fragment length (in subsampled snapshots)");
  init_treenode_info_data(list_halos_new,SID_FARG delta_fragmented,SID_DOUBLE,"delta_fragmented [Gyrs]");
  init_treenode_info_data(list_halos_new,SID_FARG Mvir_parent,     SID_DOUBLE,"M_parent [h^{-1] M_sol]");
  init_treenode_info_data(list_halos_new,SID_FARG Mvir,            SID_DOUBLE,"M_vir    [h^{-1] M_sol]");
  init_treenode_info_data(list_halos_new,SID_FARG frag_index,      SID_INT,   "Fragment's starting index");
  init_treenode_info_data(list_halos_new,SID_FARG frag_snapshot,   SID_INT,   "Fragment's starting snapshot");

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_fragmented(current_halo)){
           int flag_new=FALSE;
           if(current_halo->n_progenitors==0)
              flag_new=TRUE;
           tree_markers_info markers;
           find_treenode_markers(trees,current_halo,&markers);
           if(check_treenode_if_central(current_halo)){
              add_to_treenode_list(list_halos_central,current_halo);
              if(flag_new)
                 add_to_treenode_list(list_halos_new_central,current_halo);
           }
           else{
              add_to_treenode_list(list_halos_substructure,current_halo);
              if(flag_new)
                 add_to_treenode_list(list_halos_new_substructure,current_halo);
           }
           add_to_treenode_list(list_halos,current_halo);
           if(flag_new){
              add_to_treenode_list(list_halos_new,current_halo);
              frag_snapshot[i_list]   =trees->snap_list[current_halo->snap_tree];
              frag_index[i_list]      =current_halo->file_index;
              frag_length[i_list]     =markers.branch_root->snap_tree-current_halo->snap_tree+1;
              delta_fragmented[i_list]=fetch_treenode_delta_t(trees,markers.branch_root,current_halo)/S_PER_GYR;
              Mvir[i_list]            =fetch_treenode_Mvir   (trees,current_halo);
              Mvir_parent[i_list]     =fetch_treenode_Mvir   (trees,current_halo->parent);
              i_list++;
           }
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write data files
  write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_markers   (trees,filename_out_root,list_halos);
  write_treenode_list_properties(trees,filename_out_root,list_halos_new);
  write_treenode_list_markers   (trees,filename_out_root,list_halos_new);
  write_treenode_list_data      (trees,filename_out_root,list_halos_new);

  // Write histograms
  write_treenode_list_hist(trees,filename_out_root,list_halos);
  write_treenode_list_hist(trees,filename_out_root,list_halos_central);
  write_treenode_list_hist(trees,filename_out_root,list_halos_substructure);
  write_treenode_list_hist(trees,filename_out_root,list_halos_new);
  write_treenode_list_hist(trees,filename_out_root,list_halos_new_central);
  write_treenode_list_hist(trees,filename_out_root,list_halos_new_substructure);

  // Clean-up
  free_treenode_list(&list_halos);
  free_treenode_list(&list_halos_central);
  free_treenode_list(&list_halos_substructure);
  free_treenode_list(&list_halos_new);
  free_treenode_list(&list_halos_new_central);
  free_treenode_list(&list_halos_new_substructure);

  SID_log("Done.",SID_LOG_CLOSE);
}

