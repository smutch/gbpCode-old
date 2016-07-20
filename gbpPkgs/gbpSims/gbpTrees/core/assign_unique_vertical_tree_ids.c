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

void assign_unique_vertical_tree_ids(tree_info *trees,tree_node_info *halo){

  // Fetch this halo's properties
  halo_properties_SAGE_info *halo_properties;
  tree_node_info            *group_halo_first=NULL;
  tree_node_info            *group_halo_next =NULL;
  if(halo->parent_top==NULL)
     halo_properties=&(trees->group_properties_SAGE[halo->snap_tree][halo->neighbour_index]);
  else{
     halo_properties =&(trees->subgroup_properties_SAGE[halo->snap_tree][halo->neighbour_index]);
     group_halo_first=halo->parent_top->substructure_first;
     group_halo_next =halo->substructure_next;
  }  

  // Set the snap_num to the sampled snapshot number
  halo_properties->snap_num=halo->snap_tree;

  // Set descendant and progenitor ids
  halo_properties->descendant      =construct_unique_vertical_tree_id(halo->descendant);
  halo_properties->progenitor_first=construct_unique_vertical_tree_id(halo->progenitor_first);
  halo_properties->progenitor_next =construct_unique_vertical_tree_id(halo->progenitor_next);

  // Set group ids
  halo_properties->group_halo_first=construct_unique_vertical_tree_id(group_halo_first);
  halo_properties->group_halo_next =construct_unique_vertical_tree_id(group_halo_next);
  
}

