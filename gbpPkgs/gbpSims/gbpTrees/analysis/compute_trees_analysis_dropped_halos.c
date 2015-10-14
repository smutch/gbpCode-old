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

// Halos are processed here in such a way as to group the output by branches
void compute_trees_analysis_dropped_halos(tree_info *trees,char *filename_out_root_in,int i_type,double logM_min,double dlogM,int n_logM){

  // Compute merger rates ...
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

  SID_log("Performing dropped %sgroup analysis...",SID_LOG_OPEN|SID_LOG_TIMER,group_prefix);

  // Count the number of dropped halos and the number of branches with dropped halos
  SID_log("Counting dropped halos and branches...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_halos_dropped   =0;
  int n_branches_dropped=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=neighbour_list_start[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_DROPPED)){
           // Scan the branch ...
           // ... first, find out how many dropped instances are in the descendant line.  This
           //     is done so we can avoid double (or triple, etc) counting branches.  This will
           //     be a bit slow, since we will be rescanning halos up to n_snap times.
           int n_dropped_i=0;
           tree_node_info *current_progenitor=current_halo->descendant;
           while(current_progenitor!=NULL){
              if(current_progenitor->halo_ID!=current_halo->halo_ID) break;
              n_dropped_i+=check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_DROPPED);
              current_progenitor=current_progenitor->descendant;
           }

           // If there are no dropped instances in the descendant line, then
           //    we wont be double counting the halos in the progenitor line.
           if(n_dropped_i==0){
              n_dropped_i=0;
              current_progenitor=current_halo;
              while(current_progenitor!=NULL){
                 int flag_dropped=check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_DROPPED);
                 n_dropped_i+=flag_dropped;
                 current_progenitor=current_progenitor->progenitor_first;
                 // Because we may have reordered the progenitors, we can't trust that
                 //    the first progenitor is the one with the current halo's ID.  We 
                 //    have to scan all the progenitors as a result
                 if(current_progenitor!=NULL){
                    while(current_progenitor->halo_ID!=current_halo->halo_ID){
                       current_progenitor=current_progenitor->progenitor_next;
                       if(current_progenitor==NULL) break;
                    }
                 }
              }
              n_halos_dropped+=n_dropped_i;
              n_branches_dropped++;
           }
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  char filename_root_branches[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,   "%s_dropped_halos",   filename_out_root);
  sprintf(filename_root_branches,"%s_dropped_branches",filename_out_root);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  init_treenode_list("dropped",n_halos_dropped,&list_halos);
  int    *snap_found;   init_treenode_info_data(list_halos,SID_FARG snap_found,   SID_INT,   "Snapshot when found");
  int    *snap_lost;    init_treenode_info_data(list_halos,SID_FARG snap_lost,    SID_INT,   "Snapshot when lost");
  double *delta_normed; init_treenode_info_data(list_halos,SID_FARG delta_normed, SID_DOUBLE,"delta_dropped/t_dyn(z)");
  double *delta_dropped;init_treenode_info_data(list_halos,SID_FARG delta_dropped,SID_DOUBLE,"delta_dropped [Gyrs]");
  int    *branch_number;init_treenode_info_data(list_halos,SID_FARG branch_number,SID_INT,   "branch_number");

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_halos_dropped   =0;
  int i_branches_dropped=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=neighbour_list_start[i_snap];
     double                     z=(double)trees->z_list[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_DROPPED)){
           // Scan the branch ...
           // ... first, find out how many dropped instances are in the descendant line.  This
           //     is done so we can avoid double (or triple, etc) counting branches.  This will
           //     be a bit slow, since we will be rescanning halos up to n_snap times.
           int n_dropped_i=0;
           tree_node_info *current_progenitor=current_halo->descendant;
           while(current_progenitor!=NULL){
              if(current_progenitor->halo_ID!=current_halo->halo_ID) break;
              n_dropped_i+=check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_DROPPED);
              current_progenitor=current_progenitor->descendant;
           }

           // If there are no dropped instances in the descendant line, then
           //    we wont be double counting the halos in the progenitor line.
           if(n_dropped_i==0){
              current_progenitor=current_halo;
              while(current_progenitor!=NULL){
                 int flag_dropped=check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_DROPPED);
                 if(flag_dropped){
                    // Add an dropped halo instance to the list
                    add_to_treenode_list(list_halos,current_progenitor);
                    snap_found   [i_halos_dropped]=fetch_treenode_snapshot(trees,current_progenitor->descendant);
                    snap_lost    [i_halos_dropped]=fetch_treenode_snapshot(trees,current_progenitor);
                    delta_normed [i_halos_dropped]=fetch_treenode_delta_t (trees,current_progenitor->descendant,current_progenitor)/t_dyn_z(z,trees->cosmo);
                    delta_dropped[i_halos_dropped]=fetch_treenode_delta_t (trees,current_progenitor->descendant,current_progenitor)/S_PER_GYR;
                    branch_number[i_halos_dropped]=i_branches_dropped;
                    i_halos_dropped++;
                 }
                 current_progenitor=current_progenitor->progenitor_first;
                 // Because we may have reordered the progenitors, we can't trust that
                 //    the first progenitor is the one with the current halo's ID.  We 
                 //    have to scan all the progenitors as a result
                 if(current_progenitor!=NULL){
                    while(current_progenitor->halo_ID!=current_halo->halo_ID){
                       current_progenitor=current_progenitor->progenitor_next;
                       if(current_progenitor==NULL) break;
                    }
                 }
              }
              i_branches_dropped++;
           }
        }
        current_halo=current_halo->next_neighbour;
     }
  }

  // Finalize
  finalize_treenode_list(trees,list_halos);
  SID_log("Done.",SID_LOG_CLOSE);

  // Write the files
  write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_data      (trees,filename_out_root,list_halos);
  //write_treenode_list_hist      (trees,filename_out_root,list_halos,logM_min,dlogM,n_logM);

  // Clean-up
  free_treenode_list(&list_halos);

  SID_log("Done.",SID_LOG_CLOSE);
}

