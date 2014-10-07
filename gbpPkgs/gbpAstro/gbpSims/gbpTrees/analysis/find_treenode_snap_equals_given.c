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

int find_treenode_snap_equals_given(tree_info *trees,tree_node_info *halo,int snap_tree_given,tree_node_info **treenode_return){
   if(halo!=NULL){
      tree_node_info *treenode_current=halo;
      if(treenode_current->snap_tree>snap_tree_given){
         while(treenode_current!=NULL && treenode_current->snap_tree>snap_tree_given)
            treenode_current=treenode_current->progenitor_first;
      }
      else if(treenode_current->snap_tree<snap_tree_given){
         while(treenode_current!=NULL && treenode_current->snap_tree<snap_tree_given)
            treenode_current=treenode_current->descendant;
      }
      if(check_treenode_if_snap_equals_given(treenode_current,snap_tree_given)){
         (*treenode_return)=treenode_current;
         return(TRUE);
      }
   }
   (*treenode_return)=NULL;
   return(FALSE);
}

