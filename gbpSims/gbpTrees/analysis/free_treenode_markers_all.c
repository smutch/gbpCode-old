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

int free_treenode_markers_all(tree_info *trees,tree_markers_info ***markers){
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
      SID_free(SID_FARG (*markers)[i_snap]);
   SID_free(SID_FARG (*markers));
}

