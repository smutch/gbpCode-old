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

void assign_depth_first_index_recursive(tree_node_info *tree,int *depth_first_index){
  tree_node_info *current;

  // Set and increment the tree index
  tree->depth_first_index=(*depth_first_index)++;

  // Walk the tree
  current=tree->progenitor_first;
  while(current!=NULL){
    assign_depth_first_index_recursive(current,depth_first_index);
    current=current->progenitor_next;
  }
}

