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

void update_trees_lookup(tree_info *trees,int i_file){
   int             i_halo;
   tree_node_info *current;
   i_halo =0;
   current=trees->first_neighbour_groups[i_file];
   while(current!=NULL){
      trees->group_indices[i_file%trees->n_wrap_lookup][i_halo]=current->file_index;
      trees->group_array[i_file%trees->n_wrap_lookup][i_halo]  =current;
      i_halo++;
      current=current->next_neighbour;
   }
   i_halo =0;
   current=trees->first_neighbour_subgroups[i_file];
   while(current!=NULL){
      trees->subgroup_indices[i_file%trees->n_wrap_lookup][i_halo]=current->file_index;
      trees->subgroup_array[i_file%trees->n_wrap_lookup][i_halo]  =current;
      i_halo++;
      current=current->next_neighbour;
   }
}

