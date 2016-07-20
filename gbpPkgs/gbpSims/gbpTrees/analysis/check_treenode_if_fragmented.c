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

int check_treenode_if_fragmented(tree_node_info *halo){
   if(halo!=NULL)
      return(check_mode_for_flag(halo->tree_case,TREE_CASE_FRAGMENTED_STRAYED)||
             check_mode_for_flag(halo->tree_case,TREE_CASE_FRAGMENTED_NORMAL) ||
             check_mode_for_flag(halo->tree_case,TREE_CASE_FRAGMENTED_OTHER));
   return(FALSE);
}

