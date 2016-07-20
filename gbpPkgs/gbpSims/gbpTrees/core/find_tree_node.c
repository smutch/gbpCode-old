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

int find_tree_node(tree_info *trees,int node_file,int node_index,int group_mode,tree_node_info **found_node){
// group_mode==TRUE if we are looking for a group, FALSE for substructure
   if(node_file>=0 && node_index>=0){
      tree_node_info **halo_array;
      int             *halo_indices;
      int              index_index;
      int              n_halos;
      int              i_wrap=node_file%trees->n_wrap_lookup;
      if(group_mode==TRUE){
         n_halos     =trees->n_groups_snap_local[node_file];
         halo_indices=trees->group_indices[i_wrap];
         halo_array  =trees->group_array[i_wrap];
      }
      else{
         n_halos     =trees->n_subgroups_snap_local[node_file];
         halo_indices=trees->subgroup_indices[i_wrap];
         halo_array  =trees->subgroup_array[i_wrap];
      }
      index_index=find_index_int(halo_indices,node_index,n_halos,NULL); // halo_indices is sorted by construction
      while(halo_indices[index_index]<node_index && (index_index<(n_halos-1))) index_index++;
      if(halo_indices[index_index]!=node_index){
         (*found_node)=NULL;
         return(FALSE);
      }
      else{
         (*found_node)=halo_array[index_index];
         return(TRUE);
      }
   }
   else{
      (*found_node)=NULL;
      return(TRUE);
   }
}

