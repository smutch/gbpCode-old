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
                            tree_node_info **joined_current_parent_top){
   if(halo!=NULL){
      tree_node_info *parent_top=halo->parent_top;
      if(parent_top!=NULL){
         tree_node_info *current_halo      =NULL;
         tree_node_info *current_parent_top=NULL;
         // Find when the halo first became substructure
         if(check_treenode_if_satellite(halo))
            (*first_became_satellite)=halo;
         else
            (*first_became_satellite)=NULL;
         current_halo=halo->progenitor_first;
         while(current_halo!=NULL){
            if(check_treenode_if_satellite(current_halo))
               (*first_became_satellite)=current_halo;
            current_halo=current_halo->progenitor_first;
         }
         // Find when the halo became part of it's current top-level parent
         (*joined_current_parent_top)=halo;
         current_halo                =halo->progenitor_first;
         current_parent_top          =parent_top->progenitor_first;
         while(current_halo!=NULL && current_parent_top!=NULL){
            tree_node_info *current_halo_parent_top=current_halo->parent_top;
            if(current_parent_top==current_halo_parent_top){
               (*joined_current_parent_top)=current_halo;
               current_halo      =current_halo->progenitor_first;
               current_parent_top=current_parent_top->progenitor_first;
            }
            else{
               current_halo      =NULL;
               current_parent_top=NULL;
            }
         }
         return(TRUE);
      }
   }
   (*first_became_satellite)   =NULL;
   (*joined_current_parent_top)=NULL;
   return(FALSE);
}

