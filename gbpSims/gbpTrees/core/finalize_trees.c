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
  if(check_mode_for_flag(mode_substructure_order,TREE_SUBSTRUCTURE_ORDER_DEFAULT))
     SID_log("Assigning substructure ordering (by particle count)...",SID_LOG_OPEN|SID_LOG_TIMER);
  else
     SID_trap_error("Invalid substructure mode (%d).",ERROR_LOGIC,mode_substructure_order);
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
  if(check_mode_for_flag(mode_progenitor_order,TREE_PROGENITOR_ORDER_DELUCIA))
     SID_log("Assigning progenitor ordering (Delucia)...",SID_LOG_OPEN|SID_LOG_TIMER);
  else if(check_mode_for_flag(mode_progenitor_order,TREE_PROGENITOR_ORDER_N_PARTICLES))
     SID_log("Assigning progenitor ordering (by particle count)...",SID_LOG_OPEN|SID_LOG_TIMER);
  else
     SID_trap_error("Invalid progenitor mode (%d).",ERROR_LOGIC,mode_progenitor_order);
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

  // ... assign depth-first-indices ...
  SID_log("Assigning depth first indices...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(int i_type=0;i_type<2;i_type++){
     tree_node_info **first_in_forest;
     int             *n_halos_forest_local;
     char             group_prefix_text[5];
     switch(i_type){
        case 0:
           first_in_forest     =trees->first_in_forest_groups;
           n_halos_forest_local=trees->n_groups_forest_local;
           sprintf(group_prefix_text,"");
           break;
        case 1:
           first_in_forest     =trees->first_in_forest_subgroups;
           n_halos_forest_local=trees->n_subgroups_forest_local;
           sprintf(group_prefix_text,"sub");
           break;
     }
     for(int i_forest=0;i_forest<trees->n_forests_local;i_forest++){
        int depth_first_index=0;
        tree_node_info *current=first_in_forest[i_forest];
        while(current!=NULL){
           if(current->descendant==NULL)
              assign_depth_first_index_vertical_recursive(current,&depth_first_index);
           current=current->next_in_forest;
        }
        if(depth_first_index!=n_halos_forest_local[i_forest])
           SID_trap_error("DFI != n_halos (i.e. %d!=%d) while processing %sgroups",ERROR_LOGIC,
                          depth_first_index,n_halos_forest_local[i_forest],group_prefix_text);
     }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}

