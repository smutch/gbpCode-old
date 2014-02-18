#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void free_trees_vertical(tree_vertical_info **tree){
  tree_vertical_node_info *next_node;
  tree_vertical_node_info *last_node;
  // Free nodes
  next_node=(*tree)->root;
  while(next_node!=NULL){
    last_node=next_node;
    next_node=last_node->next;
    SID_free(SID_FARG last_node);
  }
  // Free neighbour list arrays
  SID_free(SID_FARG (*tree)->n_neighbours);
  SID_free(SID_FARG (*tree)->neighbour_halos);
  SID_free(SID_FARG (*tree)->neighbour_halo_last);
  SID_free(SID_FARG (*tree));
}

