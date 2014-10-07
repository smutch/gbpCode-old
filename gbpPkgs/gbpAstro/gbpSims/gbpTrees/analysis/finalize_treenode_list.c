#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void finalize_treenode_list(tree_info *trees,treenode_list_info *list){
   SID_Allreduce(&(list->n_list_local),&(list->n_list),1,SID_INT,SID_SUM,SID.COMM_WORLD);
}

