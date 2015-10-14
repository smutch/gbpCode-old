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

tree_node_info *fetch_treenode_backmatch(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      if(halo->parent==NULL){
         if(trees->group_backmatch_pointers!=NULL)
            return(trees->group_backmatch_pointers[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Group backmatch pointers are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_backmatch_pointers!=NULL)
            return(trees->subgroup_backmatch_pointers[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Subgroup backmatch pointers are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
   }
   else
      return(NULL);
}

