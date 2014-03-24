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

int find_treenode_last_instance(tree_info *trees,tree_node_info *halo,tree_node_info **last_instance){
   (*last_instance)=halo;
   while(check_treenode_if_main_progenitor((*last_instance)))
      (*last_instance)=(*last_instance)->descendant;
   return(TRUE);
}

