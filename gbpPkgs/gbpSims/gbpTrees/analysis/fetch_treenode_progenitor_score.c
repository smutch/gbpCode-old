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

float fetch_treenode_progenitor_score(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      if(halo->parent==NULL){
         if(trees->group_progenitor_score!=NULL)
            return(trees->group_progenitor_score[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Group progenitor scores are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_progenitor_score!=NULL)
            return(trees->subgroup_progenitor_score[halo->snap_tree][halo->neighbour_index]);
         else
            SID_trap_error("Subgroup progenitor scores are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
   }
   else
      return(-1.);
}

