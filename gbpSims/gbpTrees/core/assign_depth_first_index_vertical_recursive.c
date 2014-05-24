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

void assign_depth_first_index_vertical_recursive(tree_node_info *halo,int *depth_first_index){

  // Set and increment the depth-first index of this halo
  halo->depth_first_index=(*depth_first_index)++;

  // Walk the tree
  tree_node_info *current_progenitor=halo->progenitor_first;
  while(current_progenitor!=NULL){
    assign_depth_first_index_vertical_recursive(current_progenitor,depth_first_index);
    current_progenitor=current_progenitor->progenitor_next;
  }
}

