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

int check_treenode_if_matched_to_emerged(tree_node_info *halo){
   if(halo!=NULL)
      return(check_mode_for_flag(halo->tree_case,TREE_CASE_MATCHED_TO_EMERGED));
   return(FALSE);
}

