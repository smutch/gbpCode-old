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

int find_treenode_branch_root(tree_info *trees,tree_node_info *halo,tree_node_info **branch_root){
   (*branch_root)=halo;
   while(check_treenode_if_main_progenitor((*branch_root)))
      (*branch_root)=(*branch_root)->descendant;
   return(TRUE);
}

