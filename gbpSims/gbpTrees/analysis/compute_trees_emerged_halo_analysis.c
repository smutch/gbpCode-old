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

void compute_trees_emerged_halo_analysis(tree_info *trees,char *filename_out_root){

  // Compute merger rates ...
  SID_log("Performing emerged halo analysis...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Count the number of emerged halos and the number of branches with emerged halos
  SID_log("Counting emerged halos and branches...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_halos_emerged                     =0;
  int n_halos_emerged_central             =0;
  int n_halos_emerged_substructure        =0;
  int n_halos_emerged_nofrags             =0;
  int n_halos_emerged_nofrags_central     =0;
  int n_halos_emerged_nofrags_substructure=0;
  int n_branches_emerged=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_branch_start(trees,current_halo)){
           // Scan the branch
           int n_emerged_i                     =0;
           int n_emerged_central_i             =0;
           int n_emerged_substructure_i        =0;
           int n_emerged_nofrags_i             =0;
           int n_emerged_nofrags_central_i     =0;
           int n_emerged_nofrags_substructure_i=0;
           tree_node_info *current_progenitor=current_halo;
           while(current_progenitor!=NULL){
              // Process emerged halos in this branch's main progenitor line
              if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_EMERGED)){
                 if(!check_treenode_if_fragmented(current_progenitor)){
                    if(check_treenode_if_central(current_progenitor))
                       n_emerged_nofrags_central_i++;
                    else
                       n_emerged_nofrags_substructure_i++;
                    n_emerged_nofrags_i++;
                 }
                 if(check_treenode_if_central(current_progenitor))
                    n_emerged_central_i++;
                 else
                    n_emerged_substructure_i++;
                 n_emerged_i++;
              }
              current_progenitor=current_progenitor->progenitor_first;
           }
           if(n_emerged_i>0){
              n_halos_emerged             +=n_emerged_i;
              n_halos_emerged_central     +=n_emerged_central_i;
              n_halos_emerged_substructure+=n_emerged_substructure_i;
              n_branches_emerged++;
           }
           if(n_emerged_nofrags_i>0){
              n_halos_emerged_nofrags             +=n_emerged_nofrags_i;
              n_halos_emerged_nofrags_central     +=n_emerged_nofrags_central_i;
              n_halos_emerged_nofrags_substructure+=n_emerged_nofrags_substructure_i;
           }
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  char filename_root_branches[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,   "%s_emerged_halos",   filename_out_root);
  sprintf(filename_root_branches,"%s_emerged_branches",filename_out_root);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  treenode_list_info *list_halos_central;
  treenode_list_info *list_halos_substructure;
  treenode_list_info *list_halos_nofrags;
  treenode_list_info *list_halos_nofrags_central;
  treenode_list_info *list_halos_nofrags_substructure;
  treenode_list_info *list_branches;
  int                *branch_number;
  int                *n_emerged_list;
  int                *n_particles;
  double             *delta_normed;
  double             *delta_emerged;
  double             *Mvir;
  double             *Mvir_descendant;
  double             *Mvir_parent;
  int                *flag_is_merger;
  int                *flag_is_fraged;
  init_treenode_list("emerged_all",                 n_halos_emerged,                     &list_halos);
  init_treenode_list("emerged_centrals",            n_halos_emerged_central,             &list_halos_central);
  init_treenode_list("emerged_substructure",        n_halos_emerged_substructure,        &list_halos_substructure);
  init_treenode_list("emerged_nofrags_all",         n_halos_emerged_nofrags,             &list_halos_nofrags);
  init_treenode_list("emerged_nofrags_centrals",    n_halos_emerged_nofrags_central,     &list_halos_nofrags_central);
  init_treenode_list("emerged_nofrags_substructure",n_halos_emerged_nofrags_substructure,&list_halos_nofrags_substructure);
  init_treenode_list("emerged_branches",            n_branches_emerged,                  &list_branches);
  init_treenode_info_data(list_halos,   SID_FARG flag_is_fraged, SID_INT,   "frag'd flag");
  init_treenode_info_data(list_halos,   SID_FARG flag_is_merger, SID_INT,   "merger flag");
  init_treenode_info_data(list_halos,   SID_FARG delta_normed,   SID_DOUBLE,"delta_emerged/t_dyn(z)");
  init_treenode_info_data(list_halos,   SID_FARG delta_emerged,  SID_DOUBLE,"delta_emerged [Gyrs]");
  init_treenode_info_data(list_halos,   SID_FARG n_particles,    SID_INT,   "n_particles");
  init_treenode_info_data(list_halos,   SID_FARG Mvir_descendant,SID_DOUBLE,"M_descendant [h^{-1] M_sol]");
  init_treenode_info_data(list_halos,   SID_FARG Mvir_parent,    SID_DOUBLE,"M_parent     [h^{-1] M_sol]");
  init_treenode_info_data(list_halos,   SID_FARG Mvir,           SID_DOUBLE,"M            [h^{-1] M_sol]");
  init_treenode_info_data(list_halos,   SID_FARG branch_number,  SID_INT,   "branch_number");
  init_treenode_info_data(list_branches,SID_FARG n_emerged_list, SID_INT,   "n_emerged");

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list     =0;
  int n_emerged_i=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     double t_dyn_z=0.1*trees->t_list[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_treenode_if_branch_start(trees,current_halo)){
           // Scan the branch
           int n_emerged_j=0;
           tree_node_info *current_progenitor=current_halo;
           while(current_progenitor!=NULL){
              // Process emerged halos in this branch's main progenitor line
              if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_EMERGED)){
                 if(!check_treenode_if_fragmented(current_progenitor)){
                    if(check_treenode_if_central(current_progenitor))
                       add_to_treenode_list(list_halos_nofrags_central,current_progenitor);
                    else
                       add_to_treenode_list(list_halos_nofrags_substructure,current_progenitor);
                    add_to_treenode_list(list_halos_nofrags,current_progenitor);
                 }
                 if(check_treenode_if_central(current_progenitor))
                    add_to_treenode_list(list_halos_central,current_progenitor);
                 else
                    add_to_treenode_list(list_halos_substructure,current_progenitor);
                 add_to_treenode_list(list_halos,current_progenitor);
                 branch_number[n_emerged_i]  =i_list;
                 delta_emerged[n_emerged_i]  =fetch_treenode_delta_t(trees,current_progenitor,current_progenitor->progenitor_first)/S_PER_GYR;
                 delta_normed[n_emerged_i]   =fetch_treenode_delta_t(trees,current_progenitor,current_progenitor->progenitor_first)/t_dyn_z;
                 n_particles[n_emerged_i]    =fetch_treenode_n_particles(trees,current_progenitor);
                 Mvir[n_emerged_i]           =fetch_treenode_Mvir (trees,current_progenitor);
                 Mvir_descendant[n_emerged_i]=fetch_treenode_Mvir (trees,current_progenitor->descendant);
                 Mvir_parent[n_emerged_i]    =fetch_treenode_Mvir (trees,current_progenitor->parent);
                 flag_is_merger[n_emerged_i] =check_mode_for_flag (current_progenitor->progenitor_first->tree_case,TREE_CASE_MERGER);
                 flag_is_fraged[n_emerged_i] =check_treenode_if_fragmented(current_progenitor);
                 n_emerged_i++;
                 n_emerged_j++;
              }
              current_progenitor=current_progenitor->progenitor_first;
           }
           if(n_emerged_j>0){
              add_to_treenode_list(list_branches,current_halo);
              n_emerged_list[i_list]=n_emerged_j;
              i_list++;
           }
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write the files
  write_treenode_list_markers   (trees,filename_out_root,list_halos);
  write_treenode_list_properties(trees,filename_out_root,list_halos);
  write_treenode_list_data      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_central);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_substructure);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_nofrags);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_nofrags_central);
  write_treenode_list_hist      (trees,filename_out_root,list_halos_nofrags_substructure);
  write_treenode_list_markers   (trees,filename_out_root,list_branches);
  write_treenode_list_data      (trees,filename_out_root,list_branches);

  // Clean-up
  free_treenode_list(&list_halos);
  free_treenode_list(&list_halos_central);
  free_treenode_list(&list_halos_substructure);
  free_treenode_list(&list_halos_nofrags);
  free_treenode_list(&list_halos_nofrags_central);
  free_treenode_list(&list_halos_nofrags_substructure);
  free_treenode_list(&list_branches);

  SID_log("Done.",SID_LOG_CLOSE);
}

