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

int add_node_to_vertical_tree(tree_vertical_info  *tree,
                              int                  match_type,
                              int                  halo_id,
                              int                  group_id,
                              int                  descendant_id,
                              int                  halo_snap,
                              int                  descendant_snap,                      
                              halo_properties_SAGE_info           *properties){
  tree_vertical_node_info *new_node;
  tree_vertical_node_info *last_node;
  tree_vertical_node_info *next_node;
  tree_vertical_node_info *last_progenitor;
  tree_vertical_node_info *next_progenitor;
  tree_vertical_node_info *descendant_halo_list;
  tree_vertical_node_info *group_halo_list;

  // Create new node
  new_node=(tree_vertical_node_info *)SID_malloc(sizeof(tree_vertical_node_info));

  // Copy halo properties into new node
  memcpy(&(new_node->halo),properties,sizeof(halo_properties_SAGE_info));
  
  // Set ids and pointer defaults for the new node
  new_node->depth_first_index    =            -1; // Default; back-filled later
  new_node->n_progenitors        =             0;
  new_node->group_id             =      group_id;
  new_node->halo_id              =       halo_id;
  new_node->descendant_id        = descendant_id;
  new_node->descendant           =          NULL; // Set below
  new_node->progenitor_first     =          NULL; // Set below
  new_node->progenitor_next      =          NULL; // Set below
  new_node->progenitor_last      =          NULL; // Set below
  new_node->group_halo_first     =      new_node; // Set below if group_id>=0
  new_node->group_halo_next      =          NULL; // Set below
  new_node->neighbour_halo_next  =          NULL; // Set below
  new_node->next                 =          NULL;

  // Find the halo's descendant and set various pointers
  //   All used halos must be added to the neighbour list (done below) for this to work
  int flag_found_group=FALSE;
  if(halo_snap!=descendant_snap && descendant_snap>=0){
    descendant_halo_list=tree->neighbour_halos[descendant_snap];
    if(descendant_halo_list!=NULL && descendant_id>=0){
      next_node=descendant_halo_list;
      while(next_node!=NULL){
        if(next_node->halo_id==descendant_id){
          // Set descendant pointer and it's last progenitor
          new_node->descendant=next_node;
          new_node->descendant->n_progenitors++;
          // If this is a descendant's first progenitor, then initialize list ...
          if(new_node->descendant->progenitor_first==NULL){
            new_node->descendant->progenitor_first=new_node;
            new_node->descendant->progenitor_last =new_node;
          }
          // ... else add the new halo to the end of the list
          else{
            next_progenitor=new_node->descendant->progenitor_first;
            while(next_progenitor!=NULL){
              next_progenitor->progenitor_last=new_node;
              last_progenitor                 =next_progenitor;
              next_progenitor                 =next_progenitor->progenitor_next;
            }
            last_progenitor->progenitor_next=new_node;
          }
          break;
        }
        else
          next_node=next_node->neighbour_halo_next;
      }
      if(new_node->descendant==NULL)
        SID_trap_error("Could not find the descendant (%d;snap=%d) of a halo (%d;snap=%d)!",ERROR_LOGIC,halo_id,halo_snap,descendant_id,descendant_snap);
    }
  }
  else if(descendant_snap==halo_snap)
    SID_trap_error("This halo (%d;snap=%d) has the same snap as its descenant (%d;snap=%d)!",ERROR_LOGIC,halo_id,halo_snap,descendant_id,descendant_snap);

  // Find the new halo's group and then set the group pointers
  //   If group_id<0, then this halo is it's own group and we 
  //   can keep defaults and skip this
  group_halo_list=tree->neighbour_halos[halo_snap];
  if(group_halo_list!=NULL && group_id>=0){
    // Scan all the groups in this snapshot
    next_node=group_halo_list;
    while(next_node!=NULL){
      // If we've found a halo belonging to the same group as the new node...
      if(next_node->group_id==group_id){
        // ... then set the new node's group_halo_first pointer ...
        new_node->group_halo_first=next_node->group_halo_first;
        // ... and set the group_halo_next of the last halo that was added to the group ...
        last_node=next_node->group_halo_first;
        next_node=last_node->group_halo_next;
        while(next_node!=NULL){ // we exit the previous while automatically this way
          last_node=next_node;
          next_node=last_node->group_halo_next;
        }
        last_node->group_halo_next=new_node;
        flag_found_group=TRUE;
      }
      // ... else, keep searching
      else
        next_node=next_node->neighbour_halo_next;
    }
  }

  // Set tree neighbour pointers (needed for finding groups and descendants)
  if(tree->neighbour_halos[halo_snap]==NULL)
    tree->neighbour_halos[halo_snap]=new_node;
  else
    tree->neighbour_halo_last[halo_snap]->neighbour_halo_next=new_node;
  tree->neighbour_halo_last[halo_snap]=new_node;
  tree->n_neighbours[halo_snap]++;

  // Set the tree's root and leaf pointers (needed later for whole-tree processing)
  if(tree->root==NULL)
    tree->root=new_node;
  if(tree->last_leaf!=NULL)
    tree->last_leaf->next=new_node;
  tree->last_leaf=new_node;

  return(flag_found_group);
}

