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

void init_trees_vertical(int n_snaps,tree_vertical_info **tree){
  int i_search;
  (*tree)                     =(tree_vertical_info       *)SID_malloc(sizeof(tree_vertical_info));
  (*tree)->n_neighbours       =(int                      *)SID_malloc(sizeof(int)*n_snaps);
  (*tree)->neighbour_halos    =(tree_vertical_node_info **)SID_malloc(sizeof(tree_vertical_node_info *)*n_snaps);
  (*tree)->neighbour_halo_last=(tree_vertical_node_info **)SID_malloc(sizeof(tree_vertical_node_info *)*n_snaps);
  for(i_search=0;i_search<n_snaps;i_search++){
    (*tree)->n_neighbours[i_search]       =0;
    (*tree)->neighbour_halos[i_search]    =NULL;
    (*tree)->neighbour_halo_last[i_search]=NULL;
  }
  (*tree)->root     =NULL;
  (*tree)->last_leaf=NULL;
}

