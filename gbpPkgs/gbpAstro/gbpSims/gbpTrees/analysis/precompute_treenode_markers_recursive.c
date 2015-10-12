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

int precompute_treenode_markers_recursive(tree_info *trees,tree_markers_info **markers_array,tree_node_info *halo,tree_markers_info **markers_descendant){
   // note: - "markers_descendant" arrives as the descendat halo's markers but is used to return back a progenitor's markers
   //       - the initial call of this recursive function should be only on halos where halo->descendant==NULL
   //       - ** the values of M_peak MUST be calculated before this function is called.

   tree_markers_info *markers_halo=NULL;
   if(halo!=NULL){

      // Fetch pointer to this halo's markers
      markers_halo=&(markers_array[halo->snap_tree][halo->neighbour_index]);

      // This is the structure we want to populate:
      //   typedef struct tree_markers_info tree_markers_info;
      //   struct tree_markers_info{
      //     tree_node_info *branch_leaf;
      //     tree_node_info *branch_root;
      //     tree_node_info *descendant;
      //     tree_node_info *main_progenitor;
      //     tree_node_info *first_became_satellite;
      //     tree_node_info *joined_current_parent;
      //     tree_node_info *peak_mass;
      //     tree_node_info *half_peak_mass;
      //     tree_node_info *merger_33pc_remnant;
      //     tree_node_info *merger_33pc_host;
      //     tree_node_info *merger_33pc_merger;
      //     tree_node_info *merger_10pc_remnant;
      //     tree_node_info *merger_10pc_host;
      //     tree_node_info *merger_10pc_merger;
      //     double          M_peak;
      //   };

         
      // Set defaults (also: stuff that gets set here will be available to this halo's progenitors)
      //    Also: anything having to do with peak mass should already have been set.
      markers_halo->descendant     =halo->descendant;
      markers_halo->main_progenitor=halo->progenitor_first;
      if(halo->descendant==NULL || (*markers_descendant)==NULL) 
         markers_halo->branch_root=halo;
      else if(halo->halo_ID!=halo->descendant->halo_ID)
         markers_halo->branch_root=halo;
      else
         markers_halo->branch_root=(*markers_descendant)->branch_root;
      markers_halo->branch_leaf           =NULL;
      markers_halo->first_became_satellite=NULL;
      markers_halo->joined_current_parent =NULL;
      markers_halo->half_peak_mass        =NULL;
      markers_halo->merger_33pc_remnant   =NULL;
      markers_halo->merger_33pc_host      =NULL;
      markers_halo->merger_33pc_merger    =NULL;
      markers_halo->merger_10pc_remnant   =NULL;
      markers_halo->merger_10pc_host      =NULL;
      markers_halo->merger_10pc_merger    =NULL;

      // Inititialize some things if this is a leaf
      tree_node_info *first_progenitor=halo->progenitor_first;
      if(first_progenitor==NULL){
         if(check_treenode_if_satellite(halo)){
            markers_halo->first_became_satellite=halo;
            markers_halo->joined_current_parent =halo;
         }
         else{
            markers_halo->first_became_satellite=NULL;
            markers_halo->joined_current_parent =halo;
         }
         markers_halo->half_peak_mass        =halo;
         markers_halo->branch_leaf           =halo;
         markers_halo->merger_33pc_remnant   =NULL;
         markers_halo->merger_33pc_host      =NULL;
         markers_halo->merger_33pc_merger    =NULL;
         markers_halo->merger_10pc_remnant   =NULL;
         markers_halo->merger_10pc_host      =NULL;
         markers_halo->merger_10pc_merger    =NULL;
      }
      // Process halos that are not leaves
      else{   
         int                flag_is_a_merger       =FALSE;
         int                flag_has_a_primary     =FALSE;
         tree_node_info    *halo_main_progenitor   =NULL; // Main progenitor
         tree_node_info    *halo_secondary         =NULL; // Secondary halo of a merger
         tree_node_info    *halo_primary           =NULL; // Primary   halo of a merger
         tree_markers_info *markers_main_progenitor=NULL;
         int                n_p_peak_primary       =0;
         int                n_p_peak_secondary     =0;

         // Walk the tree and scan all the progenitors of this halo
         tree_node_info *current_progenitor=halo->progenitor_first;
         while(current_progenitor!=NULL){
            // Descend down the tree
            tree_markers_info *markers_exchange=markers_halo;
            precompute_treenode_markers_recursive(trees,markers_array,current_progenitor,&markers_exchange);
            tree_markers_info *markers_progenitor=markers_exchange;

            // Find the main progenitor (not necessarily the first progenitor,
            //    and not necessarily the primary halo in a merger)
            if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_MAIN_PROGENITOR)){
               markers_main_progenitor=markers_progenitor;
               halo_main_progenitor   =current_progenitor;
            }

            // Find the primary halo if this is a merger
            if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_MERGER_PRIMARY)){
               halo_primary      =current_progenitor;
               n_p_peak_primary  =current_progenitor->n_particles_peak;
               flag_has_a_primary=TRUE;
            }

            // Find the most massive secondary if this is a merger
            else if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_MERGER)){
               int n_p_peak_current=current_progenitor->n_particles_peak;
               if(n_p_peak_current>n_p_peak_secondary){
                  halo_secondary    =current_progenitor;
                  n_p_peak_secondary=n_p_peak_current;
               }
               flag_is_a_merger=TRUE;
            }
   
            // Move to the next progenitor
            current_progenitor=current_progenitor->progenitor_next;
         }

         // Sanity check
         if(flag_is_a_merger && ! flag_has_a_primary)
            SID_trap_error("A halo has been found to have a merger but a primary has not been marked.",ERROR_LOGIC);

         // Set formation pointers 
         find_treenode_formation(trees,halo,0.5,&(markers_halo->half_peak_mass));

         // Set leaf marker 
         markers_halo->branch_leaf=markers_main_progenitor->branch_leaf;

         // Set accretion markers
         if(markers_main_progenitor->first_became_satellite==NULL && check_treenode_if_satellite(halo))
            markers_halo->first_became_satellite=halo;
         else
            markers_halo->first_became_satellite=NULL;
         if(markers_main_progenitor->joined_current_parent==NULL  && (halo->parent)==(markers_halo->branch_root->parent))
            markers_halo->joined_current_parent=halo;
         else
            markers_halo->joined_current_parent=NULL;

         // Set merger pointers
         if(flag_is_a_merger){
            // Note: we do not check the return value of fetch_treenode_merger_info() 
            //    since errors from skips can generally be tollerated here.
            double zeta;
            fetch_treenode_merger_info(trees,&halo_secondary,&halo_primary,&zeta);
            if(halo_secondary!=NULL && halo_primary!=NULL){
               // Set new 3:1+ merger marker ...
               if(zeta>ONE_THIRD){
                  markers_halo->merger_33pc_remnant=halo;           // Time of remnant
                  markers_halo->merger_33pc_merger =halo_secondary; // Time of secondary peak mass
                  markers_halo->merger_33pc_host   =halo_primary;   // Time of secondary peak mass
               }
               // ... else if this is not a 3:1+ merger, propagate past 10:1 merger markers.
               else{
                  markers_halo->merger_33pc_remnant=markers_main_progenitor->merger_33pc_remnant;
                  markers_halo->merger_33pc_merger =markers_main_progenitor->merger_33pc_merger;
                  markers_halo->merger_33pc_host   =markers_main_progenitor->merger_33pc_host;
               }

               // Set new 10:1+ merger marker ...
               if(zeta>0.1){
                  markers_halo->merger_10pc_remnant=halo;
                  markers_halo->merger_10pc_merger =halo_secondary;
                  markers_halo->merger_10pc_host   =halo_primary;
               }
               // ... else if this is not a 10:1+ merger, propagate past 10:1 merger markers.
               else{
                  markers_halo->merger_10pc_remnant=markers_main_progenitor->merger_10pc_remnant;
                  markers_halo->merger_10pc_merger =markers_main_progenitor->merger_10pc_merger;
                  markers_halo->merger_10pc_host   =markers_main_progenitor->merger_10pc_host;
               }
            }
         }
         // ... else if this is not a merger, propagate past merger markers.
         else{
            markers_halo->merger_33pc_remnant=markers_main_progenitor->merger_33pc_remnant;
            markers_halo->merger_33pc_merger =markers_main_progenitor->merger_33pc_merger;
            markers_halo->merger_33pc_host   =markers_main_progenitor->merger_33pc_host;
            markers_halo->merger_10pc_remnant=markers_main_progenitor->merger_10pc_remnant;
            markers_halo->merger_10pc_merger =markers_main_progenitor->merger_10pc_merger;
            markers_halo->merger_10pc_host   =markers_main_progenitor->merger_10pc_host;
         }
      } // If halo is/is-not leaf
   }

   // Send a pointer to this halo's markers back to the calling function
   if(markers_descendant!=NULL) (*markers_descendant)=markers_halo; 
}

