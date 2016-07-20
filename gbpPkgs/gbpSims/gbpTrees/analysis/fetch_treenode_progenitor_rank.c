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

int fetch_treenode_progenitor_rank(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      tree_node_info *descendant  =halo->descendant;
      if(descendant!=NULL){
         tree_node_info *current_halo=descendant->progenitor_first;
         int             i_rank      =0;
         while(current_halo!=halo && current_halo!=NULL){
            i_rank++;
            current_halo=current_halo->progenitor_next;
         }
         if(current_halo!=halo)
            SID_trap_error("Invalid situation in fetch_treenode_progenitor_rank().",ERROR_LOGIC);
         return(i_rank);
      }
   }
   return(-1);
}

