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

double fetch_treenode_zeta(tree_info *trees,tree_node_info *secondary){
   if(secondary!=NULL){
      tree_node_info *descendant=secondary->descendant;
      if(descendant!=NULL){
         // Scan for the primary
         tree_node_info *primary=descendant->progenitor_first;
         while(primary!=NULL){
            if(check_mode_for_flag(primary->tree_case,TREE_CASE_MERGER_PRIMARY)) break;
            primary=primary->progenitor_next;
         }
         if(primary!=secondary && primary!=NULL){
            tree_markers_info *markers       =fetch_treenode_precomputed_markers(trees,secondary);
            if(markers!=NULL){
               tree_node_info *secondary_peak=markers->peak_mass;
               if(secondary_peak!=NULL){
                  tree_node_info *primary_peak=find_treenode_snap_equals_given(trees,primary,secondary_peak->snap_tree);
                  if(primary_peak!=NULL){
                     double M_p =fetch_treenode_Mvir(trees,primary_peak);
                     double M_s =fetch_treenode_Mvir(trees,secondary_peak);
                     double M_lo=MIN(M_p,M_s);
                     double M_hi=MAX(M_p,M_s);
                     return(M_lo/M_hi);
                  }
               }
            }
         }
      }
   }
   return(-1.);
}

