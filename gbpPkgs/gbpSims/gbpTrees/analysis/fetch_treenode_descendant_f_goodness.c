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

float fetch_treenode_descendant_f_goodness(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      if(halo->parent==NULL){
         if(trees->group_descendant_score!=NULL){
            float score=trees->group_descendant_score[halo->snap_tree][halo->neighbour_index];
            int   n_p  =halo->n_particles;
            return(match_score_f_goodness(score,n_p));
         }
         else
            SID_trap_error("Group descendant scores are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
      else{
         if(trees->subgroup_descendant_score!=NULL){
            float score=trees->subgroup_descendant_score[halo->snap_tree][halo->neighbour_index];
            int   n_p  =halo->n_particles;
            return(match_score_f_goodness(score,n_p));
         }
         else
            SID_trap_error("Subgroup descendant scores are not defined.  They probably have not been read.",ERROR_LOGIC);
      }
   }
   else
      return(-1.);
}

