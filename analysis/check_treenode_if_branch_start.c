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

int check_treenode_if_branch_start(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      if(halo->descendant!=NULL)
         return(halo!=halo->descendant->progenitor_first);
      else if(halo->snap_tree==(trees->n_snaps-1))
         return(TRUE);
   }
   return(FALSE);
}

