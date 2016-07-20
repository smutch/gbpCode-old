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

void precompute_treenode_markers_peak_mass_recursive(tree_info          *trees,
                                                     tree_markers_info **markers_array,
                                                     tree_node_info     *halo,
                                                     tree_node_info    **halo_peak_return,
                                                     int                *n_particles_peak_return,
                                                     double             *M_peak_return){
   if(halo!=NULL){
      // Fetch pointer to this halo's markers
      tree_markers_info *markers_halo=&(markers_array[halo->snap_tree][halo->neighbour_index]);

      // Initialize the result to zero
      markers_halo->peak_mass=NULL;
      markers_halo->M_peak   =0.;

      // Fetch needed properties of the halo
      double M_halo=fetch_treenode_M_vir(trees,halo);

      // Walk the tree
      tree_node_info *first_progenitor=halo->progenitor_first;
      if(first_progenitor!=NULL){
         // Find the progenitor line with the largest peak particle count
         tree_node_info *current_progenitor=first_progenitor;
         int             n_particles_peak  =0;
         while(current_progenitor!=NULL){
            // Descend down the tree
            tree_node_info *halo_peak_prog;
            int             n_particles_peak_prog;
            double          M_peak_prog;
            precompute_treenode_markers_peak_mass_recursive(trees,markers_array,current_progenitor,&halo_peak_prog,&n_particles_peak_prog,&M_peak_prog);
            // Propagate along the line of the progenitor that has the largest peak particle count.
            if(n_particles_peak_prog>=n_particles_peak){
               n_particles_peak       =n_particles_peak_prog;
               markers_halo->peak_mass=halo_peak_prog;
               markers_halo->M_peak   =M_peak_prog;
            }
            // Move to the next progenitor
            current_progenitor=current_progenitor->progenitor_next;
         }
         // Check if the current halo should set a new pointer
         if(halo->n_particles_peak>n_particles_peak){
            markers_halo->peak_mass=halo;
            markers_halo->M_peak   =M_halo;
         }
      }
      // If this halo is a leaf, initialize M_peak to be M_halo
      else{
         markers_halo->peak_mass=halo;
         markers_halo->M_peak   =M_halo;
      }

      // Send this halo's information back to its descendant
      if(halo_peak_return!=NULL)        (*halo_peak_return)       =markers_halo->peak_mass;
      if(n_particles_peak_return!=NULL) (*n_particles_peak_return)=halo->n_particles_peak;
      if(M_peak_return!=NULL)           (*M_peak_return)          =markers_halo->M_peak;
   }
}
