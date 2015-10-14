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

int init_precompute_treenode_markers(tree_info *trees,int mode){
   // Sanity check
   if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_GROUPS) &&
      check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS))
      SID_trap_error("Both group and subgroup mode is set in init_treenode_markers_all() when only one is needed.",ERROR_LOGIC);

   // Set-up to work with groups or subgroups
   tree_markers_info ***markers;
   int                 *n_halos_array;
   if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_GROUPS)){
      markers      =&(trees->group_markers);
      n_halos_array=trees->n_groups_snap_local;
   }
   else if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS)){
      markers      =&(trees->subgroup_markers);
      n_halos_array=trees->n_subgroups_snap_local;
   }
   else
      SID_trap_error("Neither group nor subgroup mode is set in init_treenode_markers_all() when one is needed.",ERROR_LOGIC);

   // Free the marker array if it is non-NULL
   if((*markers)!=NULL)
      free_precompute_treenode_markers(trees,mode);

   // Initialize the marker array
   (*markers)=(tree_markers_info **)SID_malloc(trees->n_snaps*sizeof(tree_markers_info *));
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
      (*markers)[i_snap]=(tree_markers_info *)SID_malloc(n_halos_array[i_snap]*sizeof(tree_markers_info));
}

