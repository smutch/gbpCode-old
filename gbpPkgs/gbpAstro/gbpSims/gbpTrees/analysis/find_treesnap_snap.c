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

int find_treesnap_snap(tree_info *trees,int snap){
   for(int i_snap_tree=0;i_snap_tree<trees->n_snaps;i_snap_tree++){
      if(trees->snap_list[i_snap_tree]==snap)
         return(i_snap_tree);
   }
   SID_trap_error("Could not find snapshot %d in the tree snap_list.",ERROR_LOGIC,snap);
}

