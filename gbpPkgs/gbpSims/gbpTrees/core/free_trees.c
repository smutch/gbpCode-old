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

void free_trees(tree_info **trees){

  SID_log("Freeing trees...",SID_LOG_OPEN);

  // Free reference trees
  if((*trees)->trees_reference!=NULL)
     free_trees(&((*trees)->trees_reference));

  // Free look-up arrays
  free_trees_lookup((*trees));

  // Free ADaPS structure
  ADaPS_free(SID_FARG (*trees)->data);

  // Free cosmology
  free_cosmo(&((*trees)->cosmo));

  // Free nodes
  int i_snap;
  for(i_snap=0;i_snap<(*trees)->n_snaps;i_snap++){
     tree_node_info *current;
     tree_node_info *next;
     if((*trees)->first_neighbour_groups!=NULL){
        next=(*trees)->first_neighbour_groups[i_snap];
        while(next!=NULL){
           current=next;
           next   =current->next_neighbour;
           SID_free(SID_FARG current);
        }
     }
     if((*trees)->first_neighbour_subgroups!=NULL){
        next=(*trees)->first_neighbour_subgroups[i_snap];
        while(next!=NULL){
           current=next;
           next   =current->next_neighbour;
           SID_free(SID_FARG current);
        }
     }
  }

  // Free match scores
  for(int i_type=0;i_type<2;i_type++){
     float **match_scores;
     if(i_type==0)
        match_scores=(*trees)->group_match_scores;
     else
        match_scores=(*trees)->subgroup_match_scores;
     if(match_scores!=NULL){
        for(i_snap=0;i_snap<(*trees)->n_snaps;i_snap++)
           SID_free(SID_FARG match_scores[i_snap]);
        SID_free(SID_FARG match_scores);
     }
  }

  // Free the various other arrays
  SID_free(SID_FARG (*trees)->snap_list);
  SID_free(SID_FARG (*trees)->a_list);
  SID_free(SID_FARG (*trees)->z_list);
  SID_free(SID_FARG (*trees)->t_list);
  SID_free(SID_FARG (*trees)->n_groups_catalog);
  SID_free(SID_FARG (*trees)->n_subgroups_catalog);
  SID_free(SID_FARG (*trees)->n_groups_snap_local);
  SID_free(SID_FARG (*trees)->n_subgroups_snap_local);
  SID_free(SID_FARG (*trees)->n_groups_forest_local);
  SID_free(SID_FARG (*trees)->n_subgroups_forest_local);
  SID_free(SID_FARG (*trees)->first_neighbour_groups);
  SID_free(SID_FARG (*trees)->first_neighbour_subgroups);
  SID_free(SID_FARG (*trees)->last_neighbour_groups);
  SID_free(SID_FARG (*trees)->last_neighbour_subgroups);
  SID_free(SID_FARG (*trees)->first_in_forest_groups);
  SID_free(SID_FARG (*trees)->first_in_forest_subgroups);
  SID_free(SID_FARG (*trees)->last_in_forest_groups);
  SID_free(SID_FARG (*trees)->last_in_forest_subgroups);
  SID_free(SID_FARG (*trees)->tree2forest_mapping_group);
  SID_free(SID_FARG (*trees)->tree2forest_mapping_subgroup);

  // Free the main structure
  SID_free(SID_FARG (*trees));

  SID_log("Done.",SID_LOG_CLOSE);
}

