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

int find_treenode_Mpeak(tree_info       *trees,
                        tree_node_info  *halo,
                        tree_node_info **halo_peak){
   (*halo_peak)=halo;
   if(halo!=NULL){
      tree_node_info *current_halo=halo;
      int             n_p_peak_search =halo->n_particles_peak;
      int             n_p_peak_current=current_halo->n_particles_peak;
      while((current_halo->progenitor_first)!=NULL && n_p_peak_search==n_p_peak_current){
         // Scan over all the progenitors.  In cases where there are more than one progenitor,
         //   the one labeled with TREE_CASE_MERGER_PRIMARY is the one that will take us back
         //   to the peak mass halo
         tree_node_info *dominant_progenitor=current_halo->progenitor_first;
         tree_node_info *current_progenitor =dominant_progenitor->progenitor_next;
         int n_progenitor=1;
         while(current_progenitor!=NULL){
            if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_MERGER_PRIMARY))
               dominant_progenitor=current_progenitor;
            n_progenitor++;
         }
         if(n_progenitor>1 && !check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_MERGER_PRIMARY))
            SID_trap_error("A halo has been found to have multiple progenitors, but none have been identified as a primary (file/index/type=%d/%d/%d).",ERROR_LOGIC,halo->snap_tree,halo->file_index,halo->tree_case);
         current_halo=dominant_progenitor;
      }
      // Follow the main progenitor line back.
      return(TRUE);
   }
   return(FALSE);
}

