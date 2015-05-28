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

void init_trees_lookup(tree_info *trees){
  if(check_mode_for_flag(trees->mode,TREE_MODE_REFERENCE))
     trees->n_wrap_lookup=trees->n_snaps;
  else
     trees->n_wrap_lookup=trees->n_wrap;
  trees->group_indices   =(int             **)SID_malloc(sizeof(int             *)*trees->n_wrap_lookup);
  trees->group_array     =(tree_node_info ***)SID_malloc(sizeof(tree_node_info **)*trees->n_wrap_lookup);
  trees->subgroup_indices=(int             **)SID_malloc(sizeof(int             *)*trees->n_wrap_lookup);
  trees->subgroup_array  =(tree_node_info ***)SID_malloc(sizeof(tree_node_info **)*trees->n_wrap_lookup);
  for(int i_wrap=0;i_wrap<trees->n_wrap_lookup;i_wrap++){
     trees->group_indices[i_wrap]   =(int             *)SID_malloc(sizeof(int)             *trees->n_groups_snap_alloc_local);
     trees->group_array[i_wrap]     =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_groups_snap_alloc_local);
     trees->subgroup_indices[i_wrap]=(int             *)SID_malloc(sizeof(int)             *trees->n_subgroups_snap_alloc_local);
     trees->subgroup_array[i_wrap]  =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_subgroups_snap_alloc_local);
  }
}

