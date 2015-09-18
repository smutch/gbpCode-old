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

// This function should only be called for halos with TREE_CASE_MERGER turned on.  It will return -1 otherwise.
double fetch_treenode_zeta(tree_info *trees,tree_node_info *halo_secondary){
   double zeta=-1;
   if(halo_secondary!=NULL){
      tree_node_info *descendant=halo_secondary->descendant;
      if(descendant!=NULL && check_mode_for_flag(halo_secondary->tree_case,TREE_CASE_MERGER)){
         tree_markers_info *markers_secondary=fetch_treenode_precomputed_markers(trees,halo_secondary);
         halo_secondary=markers_secondary->peak_mass;
         if(halo_secondary!=NULL){
            tree_node_info *halo_primary=descendant->progenitor_first;
            while(halo_primary!=NULL){
               if(check_mode_for_flag(halo_primary->tree_case,TREE_CASE_MERGER_PRIMARY)) break;
               halo_primary=halo_primary->progenitor_next;
            }
            halo_primary=find_treenode_snap_equals_given(trees,halo_primary,halo_secondary->snap_tree);
            int n_p_peak_secondary=0;
            int n_p_peak_primary  =1;
            if(halo_secondary!=NULL)
               n_p_peak_secondary=halo_secondary->n_particles_peak;
            if(halo_primary!=NULL)
               n_p_peak_primary=halo_primary->n_particles_peak;
            zeta=(double)n_p_peak_secondary/(double)n_p_peak_primary;
         }
      }
   }
   return(zeta);
}
