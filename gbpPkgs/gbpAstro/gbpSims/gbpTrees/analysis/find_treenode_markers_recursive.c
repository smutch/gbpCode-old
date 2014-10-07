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

int find_treenode_markers_recursive(tree_info *trees,tree_markers_info **markers_array,tree_node_info *halo,int flag_is_main_progenitor,double *M_return,tree_markers_info **markers_descendant){
   // note: - markers_exchange arrives as the descendat halo's markers but is used to return back a progenitor's markers
   //       - the initial call of this recursive function should be only on halos where halo->descendant==NULL

   double             M_halo      =0.;
   tree_markers_info *markers_halo=NULL;
   if(halo!=NULL){

      // Fetch pointer to this halo's markers
      markers_halo=&(markers_array[halo->snap_tree][halo->neighbour_index]);

      // Fetch needed properties of the halo
      M_halo=(trees->subgroup_properties[halo->snap_tree][halo->neighbour_index]).M_vir;

      // This is the structure we want to populate:
      //   typedef struct tree_markers_info tree_markers_info;
      //   struct tree_markers_info{
      //     int             flag_halo_is_main_progenitor;
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
      markers_halo->flag_halo_is_main_progenitor=flag_is_main_progenitor;
      markers_halo->descendant                  =halo->descendant;
      markers_halo->main_progenitor             =halo->progenitor_first;
      if(halo->descendant==NULL) 
         markers_halo->branch_root=halo;
      else
         markers_halo->branch_root=(*markers_descendant)->branch_root;
      markers_halo->branch_leaf                 =NULL;
      markers_halo->first_became_satellite      =NULL;
      markers_halo->joined_current_parent       =NULL;
      markers_halo->peak_mass                   =NULL;
      markers_halo->half_peak_mass              =NULL;
      markers_halo->merger_33pc_remnant         =NULL;
      markers_halo->merger_33pc_host            =NULL;
      markers_halo->merger_33pc_merger          =NULL;
      markers_halo->merger_10pc_remnant         =NULL;
      markers_halo->merger_10pc_host            =NULL;
      markers_halo->merger_10pc_merger          =NULL;
      markers_halo->M_peak                      =0.;

      // Inititialize some things if this is a leaf
      tree_node_info *first_progenitor=halo->progenitor_first;
      if(first_progenitor==NULL){
         if(check_treenode_if_satellite(halo))
            markers_halo->first_became_satellite=halo;
         else
            markers_halo->first_became_satellite=NULL;
         markers_halo->branch_leaf           =halo;
         markers_halo->joined_current_parent =NULL;
         markers_halo->peak_mass             =halo;
         markers_halo->half_peak_mass        =halo;
         markers_halo->merger_33pc_remnant   =NULL;
         markers_halo->merger_33pc_host      =NULL;
         markers_halo->merger_33pc_merger    =NULL;
         markers_halo->merger_10pc_remnant   =NULL;
         markers_halo->merger_10pc_host      =NULL;
         markers_halo->merger_10pc_merger    =NULL;
         markers_halo->M_peak                =M_halo;
      }
      else{   
         // Walk the tree
         int                flag_is_main_progenitor=TRUE;
         int                flag_is_a_merger       =FALSE;
         tree_node_info    *halo_MMP               =first_progenitor; // Most massive progenitor
         tree_node_info    *halo_MP                =first_progenitor; // Main progenitor
         tree_markers_info *markers_MP             =NULL;
         tree_markers_info *markers_MMP            =NULL;
         double             M_MMP                  =0.;
         double             M_MP                   =0.;
         tree_node_info *current_progenitor=halo->progenitor_first;
         while(current_progenitor!=NULL){
            // Descend down the tree
            double M_current_progenitor;
            tree_markers_info *markers_exchange=markers_halo;
            find_treenode_markers_recursive(trees,markers_array,current_progenitor,flag_is_main_progenitor,&M_current_progenitor,&markers_exchange);
            tree_markers_info *markers_progenitor=markers_exchange;

            // Find the most massive progenitor (MMP)
            if(flag_is_main_progenitor){
               markers_MP=markers_progenitor;
               halo_MP   =current_progenitor;
               M_MP      =M_current_progenitor;
            }
            else if(M_current_progenitor>M_MMP){
               markers_MMP     =markers_progenitor;
               halo_MMP        =current_progenitor;
               M_MMP           =M_current_progenitor;
               flag_is_a_merger=TRUE;
            }
   
            // Move to the next progenitor
            flag_is_main_progenitor=FALSE;
            current_progenitor=current_progenitor->progenitor_next;
         }

         // Set M_peak and related pointers
         find_treenode_formation(trees,halo,0.5,&(markers_halo->peak_mass),&(markers_halo->half_peak_mass));
         if(M_halo>(markers_MP->M_peak))
            markers_halo->M_peak=M_halo;
         else
            markers_halo->M_peak=markers_MP->M_peak;

         // Set leaf marker 
         markers_halo->branch_leaf=markers_MP->branch_leaf;

         // Set accretion markers
         if(markers_MP->first_became_satellite==NULL && check_treenode_if_satellite(halo))
            markers_halo->first_became_satellite=halo;
         else
            markers_halo->first_became_satellite=NULL;
         if(markers_MP->joined_current_parent==NULL  && (halo->parent)==(markers_halo->branch_root->parent))
            markers_halo->joined_current_parent=halo;
         else
            markers_halo->joined_current_parent=NULL;

         // Set merger pointers
         if(flag_is_a_merger){
            double M_peak_MP =markers_MP->M_peak;
            double M_peak_MMP=markers_MMP->M_peak;
            double M_ratio   =M_peak_MMP/M_peak_MP;
            // ... set last 3:1 merger marker ...
            if(M_ratio>ONE_THIRD){
               markers_halo->merger_33pc_remnant=halo;
               markers_halo->merger_33pc_host   =halo_MP;
               markers_halo->merger_33pc_merger =halo_MMP;
            }
            else{
               markers_halo->merger_33pc_remnant=markers_MP->merger_33pc_remnant;
               markers_halo->merger_33pc_host   =markers_MP->merger_33pc_host;
               markers_halo->merger_33pc_merger =markers_MP->merger_33pc_merger;
            }

            // ... set last 10:1 merger marker ...
            if(M_ratio>0.1){
               markers_halo->merger_10pc_remnant=halo;
               markers_halo->merger_10pc_host   =halo_MP;
               markers_halo->merger_10pc_merger =halo_MMP;
            }
            else{
               markers_halo->merger_10pc_remnant=markers_MP->merger_10pc_remnant;
               markers_halo->merger_10pc_host   =markers_MP->merger_10pc_host;
               markers_halo->merger_10pc_merger =markers_MP->merger_10pc_merger;
            }
         }
         else{
            markers_halo->merger_33pc_remnant=markers_MP->merger_33pc_remnant;
            markers_halo->merger_33pc_host   =markers_MP->merger_33pc_host;
            markers_halo->merger_33pc_merger =markers_MP->merger_33pc_merger;
            markers_halo->merger_10pc_remnant=markers_MP->merger_10pc_remnant;
            markers_halo->merger_10pc_host   =markers_MP->merger_10pc_host;
            markers_halo->merger_10pc_merger =markers_MP->merger_10pc_merger;
         }
      } // If halo is/is-not leaf
   }

   // Send a pointer to this halo's markers back tot he calling function
   if(M_return!=NULL)           (*M_return)          =M_halo;
   if(markers_descendant!=NULL) (*markers_descendant)=markers_halo; 

}

