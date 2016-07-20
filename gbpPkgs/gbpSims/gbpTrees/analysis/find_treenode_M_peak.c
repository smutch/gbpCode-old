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

int find_treenode_M_peak(tree_info       *trees,
                         tree_node_info  *halo,
                         tree_node_info **halo_peak){
   if(!check_mode_for_flag(trees->mode,TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK))
      SID_trap_error("Trees need to have their progenitors ordered by peak particle count for find_treenode_M_peak() to work.",ERROR_LOGIC);

   (*halo_peak)=halo;
   if(halo!=NULL){
      tree_node_info *current_halo=halo->progenitor_first;
      if(current_halo!=NULL){
         int n_p_peak_search =halo->n_particles_peak;
         int n_p_peak_current=current_halo->n_particles_peak;
         while(current_halo!=NULL && n_p_peak_search==n_p_peak_current){
            (*halo_peak)    =current_halo;
            current_halo    =current_halo->progenitor_first;
            if(current_halo!=NULL)
               n_p_peak_current=current_halo->n_particles_peak;
         }
      }
      return(TRUE);
   }
   return(FALSE);
}

