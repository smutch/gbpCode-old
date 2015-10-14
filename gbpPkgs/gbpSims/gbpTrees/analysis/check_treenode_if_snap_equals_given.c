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

int check_treenode_if_snap_equals_given(tree_node_info *halo,int snap_tree_given){
   if(halo!=NULL)
      return(halo->snap_tree==snap_tree_given);
   return(FALSE);
}

