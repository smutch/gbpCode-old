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
      int n_p_target=(int)(fraction*(double)(halo->n_particles_peak));
      tree_node_info *current_halo=halo;
      while(current_halo!=NULL){
         if(current_halo->n_particles_peak>=n_p_target)
            (*halo_formation)=current_halo;
         else break;
         // Because progenitors could be ordered many ways, find the most massive one
         tree_node_info *most_massive_progenitor=current_halo->progenitor_first;
         if(most_massive_progenitor!=NULL){
            int n_p_mm=most_massive_progenitor->n_particles_peak;
            current_halo=most_massive_progenitor->progenitor_next;
            while(current_halo!=NULL){
               int n_p_i=current_halo->n_particles_peak;
               if(n_p_i>n_p_mm){
                  most_massive_progenitor=current_halo;
                  n_p_mm=n_p_i;
               }
               current_halo=current_halo->progenitor_next;
            }
         }
         current_halo=most_massive_progenitor;
      }
      return(TRUE);
   }
   return(FALSE);
}

