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

int find_treenode_leaf(tree_info *trees,tree_node_info *halo,tree_node_info **leaf){
   (*leaf)=halo;
   tree_node_info *next=(*leaf)->progenitor_first;
   while(next!=NULL){
      (*leaf)=next;
      next   =(*leaf)->progenitor_first;
   }
   return(TRUE);
}

