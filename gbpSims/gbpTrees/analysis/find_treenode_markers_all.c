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

int find_treenode_markers_all(tree_info *trees,tree_markers_info ***markers){

   // We're gonna generate subgroup markers only
   int             *n_halos_local  =trees->n_subgroups_snap_local;
   tree_node_info **first_neighbour=trees->first_neighbour_subgroups;

   // Initialize the marker array
   if((*markers)==NULL){
      (*markers)=(tree_markers_info **)SID_malloc(trees->n_snaps*sizeof(tree_markers_info *));
      for(int i_snap=0;i_snap<trees->n_snaps;i_snap++)
         (*markers)[i_snap]=(tree_markers_info *)SID_malloc(n_halos_local[i_snap]*sizeof(tree_markers_info));
   }

   // Generate the markers starting recursively from each tree root
   SID_log("Generating markers...(%zd byte structue size)...",SID_LOG_OPEN|SID_LOG_TIMER,sizeof(tree_markers_info));
   for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
      tree_node_info *halo_current=first_neighbour[i_snap];
      while(halo_current!=NULL){
         if(halo_current->descendant==NULL)
            find_treenode_markers_recursive(trees,(*markers),halo_current,TRUE,NULL,NULL);
         halo_current=halo_current->next_neighbour;
      }
   }
   SID_Barrier(SID.COMM_WORLD);
   SID_log("Done.",SID_LOG_CLOSE);

}

