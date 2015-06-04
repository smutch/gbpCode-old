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

void compute_trees_dropped_halo_analysis(tree_info *trees,char *filename_out_root_in,int i_type,double logM_min,double dlogM,int n_logM){

  // Compute merger rates ...
  char group_prefix[8];
  if(i_type==0)
     sprintf(group_prefix,"");
  else
     sprintf(group_prefix,"sub");

  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s_%sgroup",filename_out_root_in,group_prefix);

  SID_log("Performing dropped %sgroup analysis...",SID_LOG_OPEN|SID_LOG_TIMER,group_prefix);

  // Count the number of dropped halos and the number of branches with dropped halos
  SID_log("Counting dropped halos...",SID_LOG_OPEN|SID_LOG_TIMER);
  int n_dropped             =0;
  int n_dropped_central     =0;
  int n_dropped_substructure=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;
     if(i_type==0)
        current_halo=trees->first_neighbour_groups[i_snap];
     else
        current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_DROPPED)){
           if(check_treenode_if_central(current_halo))
              n_dropped_central++;
           else if(check_treenode_if_satellite(current_halo))
              n_dropped_substructure++;
           n_dropped++;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Set filename roots
  char filename_root_halos[MAX_FILENAME_LENGTH];
  char filename_root_branches[MAX_FILENAME_LENGTH];
  sprintf(filename_root_halos,   "%s_dropped_%sgroup",         filename_out_root,group_prefix);
  sprintf(filename_root_branches,"%s_dropped_%sgroup_branches",filename_out_root,group_prefix);

  // Initialize the treenode list structures
  treenode_list_info *list_halos;
  treenode_list_info *list_halos_central;
  treenode_list_info *list_halos_substructure;
  int                *branch_number;
  int                *n_p_dropped;
  double             *delta_dropped;
  double             *delta_normed;
  double             *Mvir_parent;
  init_treenode_list("dropped_all",n_dropped,&list_halos);
  init_treenode_info_data(list_halos,SID_FARG Mvir_parent,  SID_DOUBLE,"M_parent [h^{-1] M_sol]");
  init_treenode_info_data(list_halos,SID_FARG delta_dropped,SID_DOUBLE,"delta_dropped [Gyrs]");
  init_treenode_info_data(list_halos,SID_FARG delta_normed, SID_DOUBLE,"delta_dropped/t_dyn(z)");
  init_treenode_info_data(list_halos,SID_FARG n_p_dropped,  SID_INT,   "n_particles");
  if(i_type==1){
     init_treenode_list("dropped_centrals",    n_dropped_central,     &list_halos_central);
     init_treenode_list("dropped_substructure",n_dropped_substructure,&list_halos_substructure);
  }

  // Create the list
  SID_log("Creating lists...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_list=0;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;
     double z=(double)trees->z_list[i_snap];
     if(i_type==0)
        current_halo=trees->first_neighbour_groups[i_snap];
     else
        current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        // Process each new branch
        if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_DROPPED)){
           if(check_treenode_if_central(current_halo))
              add_to_treenode_list(list_halos_central,current_halo);
           else if(check_treenode_if_satellite(current_halo))
              add_to_treenode_list(list_halos_substructure,current_halo);
           add_to_treenode_list(list_halos,current_halo);
           n_p_dropped[i_list]  =fetch_treenode_n_particles(trees,current_halo);
           Mvir_parent[i_list]  =fetch_treenode_Mvir       (trees,current_halo->parent);
           delta_dropped[i_list]=fetch_treenode_delta_t    (trees,current_halo->descendant,current_halo)/S_PER_GYR;
           delta_normed[i_list] =fetch_treenode_delta_t    (trees,current_halo->descendant,current_halo)/t_dyn_z(z,trees->cosmo);
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
  write_treenode_list_hist      (trees,filename_out_root,list_halos,logM_min,dlogM,n_logM);
  if(i_type==1){
     write_treenode_list_hist(trees,filename_out_root,list_halos_central,     logM_min,dlogM,n_logM);
     write_treenode_list_hist(trees,filename_out_root,list_halos_substructure,logM_min,dlogM,n_logM);
  }

  // Clean-up
  free_treenode_list(&list_halos);
  if(i_type==1){
     free_treenode_list(&list_halos_central);
     free_treenode_list(&list_halos_substructure);
  }
  SID_log("Done.",SID_LOG_CLOSE);
}

