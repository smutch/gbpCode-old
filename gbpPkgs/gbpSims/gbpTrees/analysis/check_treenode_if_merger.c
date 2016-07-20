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

int check_treenode_if_merger(tree_node_info *halo){
   int r_val=FALSE;
   if(halo!=NULL){
      tree_node_info *descendant=halo->descendant;
      if(descendant!=NULL){
         tree_node_info *progenitor_primary=descendant->progenitor_primary;
         r_val=(progenitor_primary!=NULL && progenitor_primary!=halo);
      }
   }
   return(r_val);
}
