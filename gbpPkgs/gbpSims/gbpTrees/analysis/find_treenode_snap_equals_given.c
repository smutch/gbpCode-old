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

int find_treenode_snap_equals_given(tree_info *trees,tree_node_info *halo,int snap_tree_given,tree_node_info **treenode_return,int progenitor_mode){
   int flag_success  =FALSE;
   (*treenode_return)=NULL;
   if(halo!=NULL){
      (*treenode_return)=halo;
      if((*treenode_return)->snap_tree>snap_tree_given){
         // Descend down the first progenitor line if the requested progenitor_mode 
         //    is the one used to store the progenitor order.  We have to do the check on 
         //    progenitor mode this way since the tree flag may have flags
         //    other than the progenitor order stored in it.
         if(((trees->mode)&(progenitor_mode))==progenitor_mode){
            while((*treenode_return)!=NULL && (*treenode_return)->snap_tree>snap_tree_given)
               (*treenode_return)=(*treenode_return)->progenitor_first;
         }
         // ... else, force a descent down the progenitor line with the largest peak particle count
         else if(progenitor_mode==TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK){
            while((*treenode_return)!=NULL && (*treenode_return)->snap_tree>snap_tree_given){
               // Look for the largest progenitor 
               tree_node_info *current_progenitor=(*treenode_return)->progenitor_first;
               tree_node_info *largest_progenitor=current_progenitor;
               if(current_progenitor!=NULL) current_progenitor=current_progenitor->progenitor_next;
               while(current_progenitor!=NULL){
                  if(current_progenitor->n_particles_peak>largest_progenitor->n_particles_peak)
                     largest_progenitor=current_progenitor;
                  current_progenitor=current_progenitor->progenitor_next;
               }
               (*treenode_return)=largest_progenitor;
            }
         }
         // ... else, force a descent down the progenitor line with the largest peak inclusive particle count
         else if(progenitor_mode==TREE_PROGENITOR_ORDER_N_PARTICLES_INCLUSIVE_PEAK){
            while((*treenode_return)!=NULL && (*treenode_return)->snap_tree>snap_tree_given){
               // Look for the largest progenitor
               tree_node_info *current_progenitor=(*treenode_return)->progenitor_first;
               tree_node_info *largest_progenitor=current_progenitor;
               if(current_progenitor!=NULL) current_progenitor=current_progenitor->progenitor_next;
               while(current_progenitor!=NULL){
                  if(current_progenitor->n_particles_inclusive_peak>largest_progenitor->n_particles_inclusive_peak)
                     largest_progenitor=current_progenitor;
                  current_progenitor=current_progenitor->progenitor_next;
               }
               (*treenode_return)=largest_progenitor;
            }
         }
         else
            SID_trap_error("Unsupported progenitor mode (%d) in find_treenode_snap_equals_given() for trees with mode (%d).",ERROR_LOGIC,
                           progenitor_mode,trees->mode);
      }
      else if((*treenode_return)->snap_tree<snap_tree_given){
         while((*treenode_return)!=NULL && (*treenode_return)->snap_tree<snap_tree_given)
            (*treenode_return)=(*treenode_return)->descendant;
      }
      if(check_treenode_if_snap_equals_given((*treenode_return),snap_tree_given))
         flag_success=TRUE;
   }
   return(flag_success);
}

