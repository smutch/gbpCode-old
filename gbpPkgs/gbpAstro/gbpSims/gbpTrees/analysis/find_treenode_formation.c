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

int find_treenode_formation(tree_info       *trees,
                            tree_node_info  *halo,
                            double           fraction,
                            tree_node_info **halo_formation){
   (*halo_formation)=halo;
   if(halo!=NULL){
      int n_p_peak_halo=halo->n_particles_peak;
      int n_p_target   =(int)(fraction*(double)n_p_peak_halo);
      tree_node_info *current_halo=halo;
      while(current_halo!=NULL){
         if(current_halo->n_particles_peak>n_p_target)
            (*halo_formation)=current_halo;
         else break;
         current_halo=current_halo->progenitor_first;
      }
      return(TRUE);
   }
   return(FALSE);
}

