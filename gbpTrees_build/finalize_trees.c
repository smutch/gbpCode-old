#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void finalize_trees(tree_info *trees,int mode){

  SID_log("Finalizing...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Parse mode
  int mode_substructure_order=TREE_SUBSTRUCTURE_ORDER_DEFAULT;
  int mode_progenitor_order  =mode;

  // ... correct substructure ordering ...
  SID_log("Assigning substructure ordering...",SID_LOG_OPEN|SID_LOG_TIMER);
  int i_snap;
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_group;
     current_group=trees->first_neighbour_groups[i_snap];
     while(current_group!=NULL){
        compute_substructure_order_recursive(current_group,NULL,mode_substructure_order);
        current_group=current_group->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // ... correct progenitor ordering ...
  SID_log("Assigning progenitor ordering...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     // ... groups ...
     tree_node_info *current_group;
     current_group=trees->first_neighbour_groups[i_snap];
     while(current_group!=NULL){
        if(current_group->descendant==NULL)
           compute_progenitor_order_recursive(current_group,NULL,mode_progenitor_order);
        current_group=current_group->next_neighbour;
     }
     // ... subgroups ...
     tree_node_info *current_subgroup;
     current_subgroup=trees->first_neighbour_subgroups[i_snap];
     while(current_subgroup!=NULL){
        if(current_subgroup->descendant==NULL)
           compute_progenitor_order_recursive(current_subgroup,NULL,mode_progenitor_order);
        current_subgroup=current_subgroup->next_neighbour;
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}

