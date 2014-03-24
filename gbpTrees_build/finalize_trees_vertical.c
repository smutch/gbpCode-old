#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void finalize_trees_vertical(tree_info *trees){

  SID_log("Finalizing...",SID_LOG_OPEN|SID_LOG_TIMER);

  // ... correct progenitor ordering ...
  SID_log("Correcting central halo masses...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        tree_node_info *current_parent =current_halo->parent;
        tree_node_info *current_central=current_parent->substructure_first;
        if(current_halo==current_central){
           halo_properties_SAGE_info *current_halo_properties  =&(trees->subgroup_properties_SAGE[current_halo->snap_tree][current_halo->neighbour_index]);
           halo_properties_SAGE_info *current_parent_properties=&(trees->group_properties_SAGE[current_parent->snap_tree][current_parent->neighbour_index]);
           current_halo_properties->M_vir=current_parent_properties->M_vir;
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // ... assign ids ...
  SID_log("Assigning IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;
     current_halo=trees->first_neighbour_groups[i_snap];
     while(current_halo!=NULL){
        assign_unique_vertical_tree_ids(trees,current_halo);
        current_halo=current_halo->next_neighbour;
     }
     current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        assign_unique_vertical_tree_ids(trees,current_halo);
        current_halo=current_halo->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}

