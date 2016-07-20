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

int fetch_treenode_substructure_rank(tree_info *trees_in,tree_node_info *halo_in){
   int rank=0; // default value for top-level substructure
   if(halo_in!=NULL){
      tree_node_info *parent=halo_in->parent;
      while(parent!=NULL){
         rank++;
         parent=parent->parent;
      }
      return(rank);
   }
   return(-1);
}

