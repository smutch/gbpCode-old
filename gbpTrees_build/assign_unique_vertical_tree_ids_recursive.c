#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void assign_unique_vertical_tree_ids_recursive(tree_vertical_node_info *tree_node,int i_tree){
  tree_vertical_node_info *current;
  
  // Set descendant and progenitor ids
  tree_node->halo.descendant      =construct_unique_vertical_tree_id(tree_node->descendant,      i_tree);
  tree_node->halo.progenitor_first=construct_unique_vertical_tree_id(tree_node->progenitor_first,i_tree);
  tree_node->halo.progenitor_next =construct_unique_vertical_tree_id(tree_node->progenitor_next, i_tree);

  // Set group ids; needs to be modified for MPI
  tree_node->halo.group_halo_first=construct_unique_vertical_tree_id(tree_node->group_halo_first,i_tree);
  tree_node->halo.group_halo_next =construct_unique_vertical_tree_id(tree_node->group_halo_next, i_tree);
  
  // Walk the tree
  current=tree_node->progenitor_first;
  while(current!=NULL){
    assign_unique_vertical_tree_ids_recursive(current,i_tree);
    current=current->progenitor_next;
  }
}

