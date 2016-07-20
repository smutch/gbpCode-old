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

double fetch_treenode_match_score(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      float match_score;
      if(halo->parent_top==NULL){
         if(trees->group_match_scores!=NULL)
            match_score=trees->group_match_scores[halo->snap_tree][halo->neighbour_index];
         else
            SID_trap_error("Group match scores are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_match_scores!=NULL)
            match_score=trees->subgroup_match_scores[halo->snap_tree][halo->neighbour_index];
         else
            SID_trap_error("Subgroup match scores are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      return((double)match_score);
   }
   return(-1.);
}

