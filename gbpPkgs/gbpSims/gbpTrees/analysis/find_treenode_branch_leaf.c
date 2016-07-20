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

int find_treenode_branch_leaf(tree_info *trees,tree_node_info *halo,tree_node_info **branch_leaf){
   if(halo!=NULL){
      (*branch_leaf)=halo;
      tree_node_info *next=(*branch_leaf)->progenitor_first;
      while(next!=NULL){
         (*branch_leaf)=next;
         next          =(*branch_leaf)->progenitor_first;
      }
      return(TRUE);
   }
   (*branch_leaf)=NULL;
   return(FALSE);
}

