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

int init_treenode_markers_all(tree_info *trees,tree_markers_info ***markers,int mode){
   // Sanity check
   if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_GROUPS) &&
      check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS))
      SID_trap_error("Both group and subgroup mode is set in init_treenode_markers_all() when only one is needed.",ERROR_LOGIC);
   // Free the marker array if it is non-NULL
   if((*markers)!=NULL)
      free_treenode_markers_all(trees,markers);
   // Initialize the marker array
   (*markers)=(tree_markers_info **)SID_malloc(trees->n_snaps*sizeof(tree_markers_info *));
   if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_GROUPS)){
      for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
         (*markers)[i_snap]=(tree_markers_info *)SID_malloc(trees->n_groups_snap_local[i_snap]*sizeof(tree_markers_info));
   }
   else if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS)){
      for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
         (*markers)[i_snap]=(tree_markers_info *)SID_malloc(trees->n_subgroups_snap_local[i_snap]*sizeof(tree_markers_info));
   }
   else
      SID_trap_error("Neither group nor subgroup mode is set in init_treenode_markers_all() when one is needed.",ERROR_LOGIC);
}

