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

int add_node_to_trees(tree_info        *trees,            // The tree datastructure
                      int               i_forest,         // Local forest index
                      int               tree_case,        // Halo's TREE_CASE BWS
                      int               n_particles,      // Number of particles in the halo
                      int               halo_snap,        // Halo's tree snapshot number
                      int               halo_index,       // Halo's file index
                      int               descendant_snap,  // Descendant's snap
                      int               descendant_index, // Descendant's snap
                      tree_node_info   *group_node,       // Pointer to the new node's group
                      tree_node_info  **new_node){        // Pointer to the new node
  int rval=TRUE;

  // Create new node
  (*new_node)=(tree_node_info *)SID_malloc(sizeof(tree_node_info));

  // Set ids and pointer defaults for the new node
  (*new_node)->n_progenitors      =           0;
  (*new_node)->n_substructures    =           0;
  (*new_node)->n_particles        = n_particles;
  (*new_node)->tree_case          =   tree_case;
  (*new_node)->snap_tree          =   halo_snap;
  (*new_node)->file_index         =  halo_index;
  if(group_node==NULL)
     (*new_node)->neighbour_index = trees->n_groups_snap_local[halo_snap];
  else
     (*new_node)->neighbour_index = trees->n_subgroups_snap_local[halo_snap];
  (*new_node)->parent             =  group_node;
  (*new_node)->substructure_first =        NULL; 
  (*new_node)->substructure_last  =        NULL; 
  (*new_node)->substructure_next  =        NULL; 
  (*new_node)->descendant         =        NULL;
  (*new_node)->progenitor_first   =        NULL;
  (*new_node)->progenitor_last    =        NULL;
  (*new_node)->progenitor_next    =        NULL;
  (*new_node)->next_neighbour     =        NULL;
  (*new_node)->next_in_forest     =        NULL;

  // Are we processing a group?
  int flag_processing_group=(group_node==NULL);

  // Process descendants
  if(descendant_snap>=0){
     if(descendant_snap==halo_snap && descendant_snap>=0){
        if(flag_processing_group)
           SID_trap_error("This group halo (snap=%d;idx=%d) has the same snap as its descenant (snap=%d;idx=%d)!",ERROR_LOGIC,
                          halo_snap,halo_index,descendant_snap,descendant_index);
        else
           SID_trap_error("This subgroup halo (snap=%d;idx=%d) has the same snap as its descenant (snap=%d;idx=%d)!",ERROR_LOGIC,
                          halo_snap,halo_index,descendant_snap,descendant_index);
     }
     if(descendant_index>=0){
        int             index_index;
        tree_node_info *descendant;
        // Find the descendant
        if(!find_tree_node(trees,descendant_snap,descendant_index,flag_processing_group,&descendant)){
          if(flag_processing_group)
             SID_trap_error("Could not find descendant group (snap=%d->%d, index=%d->%d)",ERROR_LOGIC,
                            halo_snap,descendant_snap,halo_index,descendant_index);
          else
             SID_trap_error("Could not find descendant subgroup (snap=%d->%d, index=%d->%d)",ERROR_LOGIC,
                            halo_snap,descendant_snap,halo_index,descendant_index);
        }
        if(descendant==NULL){
           if(flag_processing_group)
              SID_trap_error("A group (snap=%d;idx=%d) has been found to have an undefined descendant (snap=%d;idx=%d)!",ERROR_LOGIC,
                             halo_snap,halo_index,descendant_snap,descendant_index);
           else
              SID_trap_error("A subgroup (snap=%d;idx=%d) has been found to have an undefined descendant (snap=%d;idx=%d)!",ERROR_LOGIC,
                             halo_snap,halo_index,descendant_snap,descendant_index);
        }
        (*new_node)->descendant=descendant;
        // Update its info, pointers, etc
        descendant->n_progenitors++;
        if(descendant->progenitor_first==NULL)
          descendant->progenitor_first=(*new_node);
        else
          descendant->progenitor_last->progenitor_next=(*new_node);
        descendant->progenitor_last=(*new_node);
     }
  }

  // Set substructure pointers
  if(group_node!=NULL){
     group_node->n_substructures++;
     if(group_node->substructure_first==NULL)
        group_node->substructure_first=(*new_node);
     else
        group_node->substructure_last->substructure_next=(*new_node);
     group_node->substructure_last=(*new_node);
  }

  // Set neighbour pointers
  if(group_node==NULL){
    trees->n_groups_snap_local[halo_snap]++;
    if(trees->first_neighbour_groups[halo_snap]==NULL)
      trees->first_neighbour_groups[halo_snap]=(*new_node);
    else
      trees->last_neighbour_groups[halo_snap]->next_neighbour=(*new_node);
    trees->last_neighbour_groups[halo_snap]=(*new_node);
  }
  else{
    trees->n_subgroups_snap_local[halo_snap]++;
    if(trees->first_neighbour_subgroups[halo_snap]==NULL)
      trees->first_neighbour_subgroups[halo_snap]=(*new_node);
    else
      trees->last_neighbour_subgroups[halo_snap]->next_neighbour=(*new_node);
    trees->last_neighbour_subgroups[halo_snap]=(*new_node);
  }

  // Set forest pointers
  if(group_node==NULL){
    trees->n_groups_forest_local[i_forest]++;
    if(trees->first_in_forest_groups[i_forest]==NULL)
      trees->first_in_forest_groups[i_forest]=(*new_node);
    else
      trees->last_in_forest_groups[i_forest]->next_in_forest=(*new_node);
    trees->last_in_forest_groups[i_forest]=(*new_node);
  }
  else{
    trees->n_subgroups_forest_local[i_forest]++;
    if(trees->first_in_forest_subgroups[i_forest]==NULL)
      trees->first_in_forest_subgroups[i_forest]=(*new_node);
    else
      trees->last_in_forest_subgroups[i_forest]->next_in_forest=(*new_node);
    trees->last_in_forest_subgroups[i_forest]=(*new_node);
  }  

  return(rval);
}

