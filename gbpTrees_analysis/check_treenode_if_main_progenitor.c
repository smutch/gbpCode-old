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

int check_treenode_if_main_progenitor(tree_node_info *halo){
   tree_node_info *descendant=halo->descendant;
   if(descendant!=NULL){
      tree_node_info *main_progenitor=descendant->progenitor_first;
      return(halo==main_progenitor);
   }
   else
      return(FALSE);
}

