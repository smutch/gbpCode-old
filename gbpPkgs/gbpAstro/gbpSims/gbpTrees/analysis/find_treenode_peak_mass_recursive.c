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

int find_treenode_peak_mass_recursive(tree_info *trees,tree_markers_info **markers_array,tree_node_info *halo,int *halo_type_return,double *M_peak_return){
   double             M_halo      =0.;
   tree_markers_info *markers_halo=NULL;
   if(halo!=NULL){

      // Fetch pointer to this halo's markers
      markers_halo=&(markers_array[halo->snap_tree][halo->neighbour_index]);

      // Initialize the result to zero
      markers_halo->peak_mass=NULL;
      markers_halo->M_peak   =0.;

      // Fetch needed properties of the halo
      M_halo=fetch_treenode_Mvir(trees,halo);

      tree_node_info *first_progenitor=halo->progenitor_first;
      if(first_progenitor!=NULL){   
         // Walk the tree
         tree_node_info *primary_halo;
         int             progenitor_type;
         double          M_peak_progenitor;
         double          M_peak_primary;
         tree_node_info *current_progenitor=halo->progenitor_first;
         while(current_progenitor!=NULL){
            // Descend down the tree
            find_treenode_peak_mass_recursive(trees,markers_array,current_progenitor,&progenitor_type,&M_peak_progenitor);

            // Propagate the dominant halo's peak mass
            if(check_mode_for_flag(progenitor_type,TREE_CASE_MERGER_PRIMARY)){
               primary_halo  =current_progenitor;
               M_peak_primary=M_peak_progenitor;
            }

            // Move to the next progenitor
            current_progenitor=current_progenitor->progenitor_next;
         }
         // Perform propagation
         markers_halo->peak_mass=primary_halo;
         markers_halo->M_peak   =M_peak_progenitor;
      }
      // If this halo is a leaf, initialize M_peak to be M_halo
      else
         markers_halo->M_peak=M_halo;

      // Send this halo's information back to its descendant
      if(halo_type_return!=NULL) (*halo_type_return)=halo->tree_case;
      if(M_peak_return!=NULL)    (*M_peak_return)   =markers_halo->M_peak;
   }
}

