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
   int flag_success  =FALSE;
   (*treenode_return)=NULL;
   if(halo!=NULL){
      (*treenode_return)=halo;
      if((*treenode_return)->snap_tree>snap_tree_given){
         while((*treenode_return)!=NULL && (*treenode_return)->snap_tree>snap_tree_given)
            (*treenode_return)=(*treenode_return)->progenitor_first;
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

