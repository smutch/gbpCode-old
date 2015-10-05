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

void compute_trees_analysis_fragmented_halos(tree_info *trees,char *filename_out_root_in,int i_type,double logM_min,double dlogM,int n_logM){

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

  SID_log("Performing fragmented %sgroup analysis...",SID_LOG_OPEN|SID_LOG_TIMER,group_prefix);

  // Count the number of fragmented halos 
  SID_log("Counting fragmented halos...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_fragmented                 =0;
  int n_fragmented_new             =0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=neighbour_list_start[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_fragmented(current_halo)){
           if(current_halo->n_progenitors==0)
              n_fragmented_new++;
           n_fragmented++;
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
  treenode_list_info *list_halos_new;
  init_treenode_list("fragmented_new",n_fragmented_new,&list_halos_new);
  init_treenode_list("fragmented",    n_fragmented,    &list_halos);
  int    *frag_length;     init_treenode_info_data(list_halos_new,SID_FARG frag_length,     SID_INT,   "Fragment length (snapshots)");
  double *delta_frag_norm; init_treenode_info_data(list_halos_new,SID_FARG delta_frag_norm, SID_DOUBLE,"delta_fragmented/t_dyn(z)");
  double *delta_fragmented;init_treenode_info_data(list_halos_new,SID_FARG delta_fragmented,SID_DOUBLE,"delta_fragmented [Gyrs]");
  int    *stop_index;      init_treenode_info_data(list_halos_new,SID_FARG stop_index,      SID_INT,   "Fragment's stopping index");
  int    *stop_snapshot;   init_treenode_info_data(list_halos_new,SID_FARG stop_snapshot,   SID_INT,   "Fragment's stopping snapshot");
  int    *start_index;     init_treenode_info_data(list_halos_new,SID_FARG start_index,     SID_INT,   "Fragment's starting index");
  int    *start_snapshot;  init_treenode_info_data(list_halos_new,SID_FARG start_snapshot,  SID_INT,   "Fragment's starting snapshot");

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list    =0;
  int i_list_all=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=neighbour_list_start[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_fragmented(current_halo)){
           tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,current_halo);
           int flag_new=FALSE;
           if(current_halo->n_progenitors==0)
              flag_new=TRUE;
           int flag_central_i  =check_treenode_if_central  (current_halo);
           int flag_satellite_i=check_treenode_if_satellite(current_halo);
           if(flag_new){
              add_to_treenode_list(list_halos_new,current_halo);
              frag_length[i_list]     =fetch_treenode_snapshot(trees,markers->branch_root)-fetch_treenode_snapshot(trees,current_halo);
              delta_frag_norm[i_list] =fetch_treenode_delta_t(trees,markers->branch_root,current_halo)/t_dyn_z(z,trees->cosmo);
              delta_fragmented[i_list]=fetch_treenode_delta_t(trees,markers->branch_root,current_halo)/S_PER_GYR;
              stop_index[i_list]      =markers->branch_root->file_index;
              stop_snapshot[i_list]   =trees->snap_list[markers->branch_root->snap_tree];
              start_index[i_list]     =current_halo->file_index;
              start_snapshot[i_list]  =trees->snap_list[current_halo->snap_tree];
              i_list++;
           }
           add_to_treenode_list(list_halos,current_halo);
           i_list_all++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }

  // Finalize
  finalize_treenode_list(trees,list_halos);
  finalize_treenode_list(trees,list_halos_new);
  SID_log("Done.",SID_LOG_CLOSE);

  // Write data files
  write_treenode_list_properties(trees,filename_out_root,list_halos_new);
  write_treenode_list_data      (trees,filename_out_root,list_halos_new);
  write_treenode_list_properties(trees,filename_out_root,list_halos);

  // Write histograms
  //write_treenode_list_hist(trees,filename_out_root,list_halos_new,logM_min,dlogM,n_logM);
  //write_treenode_list_hist(trees,filename_out_root,list_halos,    logM_min,dlogM,n_logM);

  // Clean-up
  free_treenode_list(&list_halos);
  free_treenode_list(&list_halos_new);
  SID_log("Done.",SID_LOG_CLOSE);
}

