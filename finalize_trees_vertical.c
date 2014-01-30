#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

void finalize_trees_vertical(tree_info **trees,
                             int        *n_halos_tree_local,
                             int         n_trees,
                             int         n_snaps,
                             int         progenitor_mode){
  int              i_tree,i_snap;
  tree_node_info  *current;
  tree_node_info  *next;
  int              progenitor_score;
  int              depth_first_index;

  SID_log("Finalizing...",SID_LOG_OPEN|SID_LOG_TIMER);

  // ... correct group halo ordering ...
  SID_log("Assigning group ordering...",SID_LOG_OPEN);
  for(i_tree=0;i_tree<n_trees;i_tree++){
    for(i_snap=0;i_snap<n_snaps;i_snap++)
      assign_group_subgroup_order(trees[i_tree],i_snap,progenitor_mode);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // ... correct progenitor ordering ...
  SID_log("Correcting progenitor ordering...",SID_LOG_OPEN);
  for(i_tree=0;i_tree<n_trees;i_tree++){
    current=trees[i_tree]->root;
    while(current!=NULL){
      if(current->descendant==NULL){
        progenitor_score=0;
        assign_progenitor_order_recursive(current,&progenitor_score,progenitor_mode);
      }
      current=current->next;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // ... assign depth-first-indices ...
  SID_log("Assigning depth first indices...",SID_LOG_OPEN);
  for(i_tree=0;i_tree<n_trees;i_tree++){
    depth_first_index=0;
    current=trees[i_tree]->root;
    while(current!=NULL){
      if(current->descendant==NULL)
        assign_depth_first_index_recursive(current,&depth_first_index);
      current=current->next;
    }
    if(depth_first_index!=n_halos_tree_local[i_tree])
      SID_trap_error("DFI != n_halos (i.e. %d!=%d)",ERROR_LOGIC,depth_first_index,n_halos_tree_local[i_tree]);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // ... assign ids ...
  SID_log("Assigning IDs...",SID_LOG_OPEN);
  for(i_tree=0;i_tree<n_trees;i_tree++){
    current=trees[i_tree]->root;
    while(current!=NULL){
      if(current->descendant==NULL)
        assign_unique_vertical_tree_ids_recursive(current,i_tree);
      current=current->next;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}

