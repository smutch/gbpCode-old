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

int find_treenode_accretion(tree_info       *trees,
                            tree_node_info  *halo,
                            tree_node_info **first_accretion_progenitor,
                            tree_node_info **last_accretion_progenitor){
   tree_node_info *parent=halo->parent;
   if(halo!=NULL && parent!=NULL){
      tree_node_info *current_halo;
      tree_node_info *current_parent;
      (*first_accretion_progenitor) =halo;
      (*last_accretion_progenitor)  =halo;
      // Find when the halo first became substructure
      current_halo=halo->progenitor_first;
      while(current_halo!=NULL){
         tree_node_info *current_halo_parent=current_halo->parent;
         if(current_halo_parent->substructure_first!=current_halo)
            (*first_accretion_progenitor)=current_halo;
         current_halo  =current_halo->progenitor_first;
      }
      // Find when the halo first joined it's current parent
      current_halo  =halo->progenitor_first;
      current_parent=parent->progenitor_first;
      while(current_halo!=NULL && current_parent!=NULL){
         tree_node_info *current_halo_parent=current_halo->parent;
         if(current_parent==current_halo_parent)
            (*last_accretion_progenitor)=current_halo;
         current_halo  =current_halo->progenitor_first;
         current_parent=current_parent->progenitor_first;
      }
      return(TRUE);
   }
   (*first_accretion_progenitor)=NULL;
   (*last_accretion_progenitor) =NULL;
   return(FALSE);
}

