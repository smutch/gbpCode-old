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

int precompute_treenode_markers(tree_info *trees,int mode){

   // Sanity check
   if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_GROUPS) && 
      check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS))
      SID_trap_error("Both group and subgroup mode flags are set when only one is allowed in precompute_treenode_markers().",ERROR_LOGIC);

   // Initialize the marker array
   init_precompute_treenode_markers(trees,mode);

   // First, process the reference trees if present
   if(trees->trees_reference!=NULL)
      precompute_treenode_markers(trees->trees_reference,mode);

   char                group_text[32];
   int                 flag_process_groups;
   int                *n_halos_snap_local   =NULL;
   tree_node_info    **first_neighbours=NULL;
   tree_markers_info **markers;
   if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_GROUPS)){
      sprintf(group_text,"groups");
      n_halos_snap_local =trees->n_groups_snap_local;
      first_neighbours   =trees->first_neighbour_groups;
      markers            =trees->group_markers;
      flag_process_groups=TRUE;
   }
   else if(check_mode_for_flag(mode,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS)){
      sprintf(group_text,"subgroups");
      n_halos_snap_local =trees->n_subgroups_snap_local;
      first_neighbours   =trees->first_neighbour_subgroups;
      markers            =trees->subgroup_markers;
      flag_process_groups=FALSE;
   }
   else
      SID_trap_error("Neither group nor subgroup mode flags are set when only one is allowed in precompute_treenode_markers().",ERROR_LOGIC);

   // Generate the markers starting recursively from each tree root
   SID_log("Generating %s markers...(%zd byte structue size)...",SID_LOG_OPEN|SID_LOG_TIMER,group_text,sizeof(tree_markers_info));
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
      SID_log("Processing snapshot #%03d of %03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap+1,trees->n_snaps);
      tree_node_info *halo_current=first_neighbours[i_snap];
      while(halo_current!=NULL){
         if(halo_current->descendant==NULL){
            find_treenode_peak_mass_recursive(trees,markers,halo_current,NULL,NULL,NULL);
            find_treenode_markers_recursive  (trees,markers,halo_current,TRUE,NULL,NULL);
         }
         halo_current=halo_current->next_neighbour;
      }
      SID_log("Done.",SID_LOG_CLOSE);
   }
   SID_Barrier(SID.COMM_WORLD);
   SID_log("Done.",SID_LOG_CLOSE);

}

