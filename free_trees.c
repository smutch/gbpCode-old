#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void free_trees(tree_info **trees){

  SID_log("Freeing trees...",SID_LOG_OPEN);
  int i_snap;
  for(i_snap=0;i_snap<(*trees)->n_snaps;i_snap++){
     tree_node_info *current;
     tree_node_info *next;
     next=(*trees)->first_neighbour_groups[i_snap];
     while(next!=NULL){
        current=next;
        next   =current->next_neighbour;
        SID_free(SID_FARG current);
     }
     next=(*trees)->first_neighbour_subgroups[i_snap];
     while(next!=NULL){
        current=next;
        next   =current->next_neighbour;
        SID_free(SID_FARG current);
     }
  }
  SID_free(SID_FARG (*trees)->snap_list);
  SID_free(SID_FARG (*trees)->z_list);
  SID_free(SID_FARG (*trees)->t_list);
  SID_free(SID_FARG (*trees)->n_groups_snap_local);
  SID_free(SID_FARG (*trees)->n_subgroups_snap_local);
  SID_free(SID_FARG (*trees)->first_neighbour_groups);
  SID_free(SID_FARG (*trees)->first_neighbour_subgroups);
  SID_free(SID_FARG (*trees)->last_neighbour_groups);
  SID_free(SID_FARG (*trees)->last_neighbour_subgroups);
  SID_free(SID_FARG (*trees)->first_in_forest_groups);
  SID_free(SID_FARG (*trees)->first_in_forest_subgroups);
  SID_free(SID_FARG (*trees)->last_in_forest_groups);
  SID_free(SID_FARG (*trees)->last_in_forest_subgroups);
  SID_free(SID_FARG (*trees)->n_groups_forest_local);
  SID_free(SID_FARG (*trees)->n_subgroups_forest_local);
  SID_free(SID_FARG (*trees));
  SID_log("Done.",SID_LOG_CLOSE);

}

