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
                            tree_node_info **first_became_satellite,
                            tree_node_info **joined_current_parent){
   if(halo!=NULL){
      tree_node_info *parent=halo->parent;
      if(parent!=NULL){
         tree_node_info *current_halo;
         tree_node_info *current_parent;
         // Find when the halo first became substructure
         if(check_treenode_if_satellite(halo))
            (*first_became_satellite)=halo;
         else
            (*first_became_satellite)=NULL;
         current_halo=halo->progenitor_first;
         while(current_halo!=NULL){
            tree_node_info *current_halo_parent=current_halo->parent;
            if(check_treenode_if_satellite(current_halo))
               (*first_became_satellite)=current_halo;
            current_halo  =current_halo->progenitor_first;
         }
         // Find when the halo became part of it's current parent
         (*joined_current_parent)=halo;
         current_halo  =halo->progenitor_first;
         current_parent=parent->progenitor_first;
         while(current_halo!=NULL && current_parent!=NULL){
            tree_node_info *current_halo_parent=current_halo->parent;
            if(current_parent==current_halo_parent){
               (*joined_current_parent)=current_halo;
               current_halo  =current_halo->progenitor_first;
               current_parent=current_parent->progenitor_first;
            }
            else{
               current_halo  =NULL;
               current_parent=NULL;
            }
         }
         return(TRUE);
      }
   }
   (*first_became_satellite)=NULL;
   (*joined_current_parent) =NULL;
   return(FALSE);
}

