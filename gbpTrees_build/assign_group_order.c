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

void assign_group_order(tree_info *trees,int mode){
  int              n_subgroups_snap_local;
  int              i_halo,j_halo,k_halo;
  int             *group_ids;
  int              largest_group;
  int              n_in_group;
  tree_node_info **neighbours;
  tree_node_info  *first_neighbour;
  tree_node_info  *current;
  tree_node_info  *new_first;
  tree_node_info  *new_last;
  int             *halo_score;
  size_t          *group_ids_index;
  size_t          *halo_score_index;
  size_t          *halo_score_rank;

  // Allocate a temporary score array large enough
  //   even for the largest locally stored group
  int *progenitor_score;
  progenitor_score=(int *)SID_malloc(sizeof(int)*trees->n_substructure_max_local);

  // Loop over each snapshot
  int i_snap;
  for(i_snap=0;i_snap<tree->n_snaps;i_snap++){
     // Loop over each subgroup in the snapshot
     tree_node_info *current_group;
     current_group=trees->first_neighbour_groups[i_snap];
     while(current_group!=NULL){
        assign_substructure_order_recursive(current_group,&progenitor_score,mode);
        current_group=current_group->next_neighbour;
     }
  }

  // Clean-up
  SID_free(SID_FARG progenitor_score);

}
/*





     n_subgroups_snap_local=tree->n_subgroups_snap_local[i_snap];
     if(n_subgroups_snap_local>0){
       // Initialize some temporary arrays
       group_ids      =(int             *)SID_malloc(sizeof(int)*n_subgroups_snap_local);
       neighbours     =(tree_node_info **)SID_malloc(sizeof(tree_node_info)*n_subgroups_snap_local);
       first_neighbour=tree->first_neighbour_subgroups[i_snap];

       // Sort halos by group_id
       i_halo =0;
       current=first_neighbour;
       while(current!=NULL){
         if(i_halo>=n_subgroups_snap_local)
           SID_trap_error("There's a problem with the number of neighbours in assign_group_subgroup_order (%d>=%d)!",ERROR_LOGIC,i_halo,n_subgroups_snap_local);
         neighbours[i_halo]=current;
         group_ids[i_halo] =current->group_id;
         current           =current->neighbour_halo_next;
         i_halo++;
       }
       if(i_halo!=n_subgroups_snap_local)
         SID_trap_error("There's a problem with the number of neighbours in assign_group_subgroup_order (%d!=%d)!",ERROR_LOGIC,i_halo,n_subgroups_snap_local);
       merge_sort(group_ids,(size_t)n_subgroups_snap_local,&group_ids_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

       // Find highest group score and initialize a temporary array
       largest_group=0;
       i_halo       =0;
       while(i_halo<n_subgroups_snap_local){
         n_in_group=0;
         while(group_ids[group_ids_index[i_halo+n_in_group]]==group_ids[group_ids_index[i_halo]]){ 
           n_in_group++;
           if((i_halo+n_in_group)>=n_subgroups_snap_local)
             break;
         }
         largest_group=MAX(largest_group,n_in_group);
         i_halo      +=n_in_group;
       }

       // Scan neighbour list (now ordered by ID) and order by score within each group
       if(largest_group>1){
         halo_score=(int *)SID_malloc(sizeof(int)*largest_group);
         i_halo   =0;
         while(i_halo<n_subgroups_snap_local){
           // Create a list of halo scores for each group
           n_in_group=0;
           while(group_ids[group_ids_index[i_halo+n_in_group]]==group_ids[group_ids_index[i_halo]]){
             halo_score[n_in_group]=0;
             compute_progenitor_score_recursive(neighbours[group_ids_index[i_halo+n_in_group]],&(halo_score[n_in_group]),mode);
             n_in_group++;
             if((i_halo+n_in_group)>=n_subgroups_snap_local)
               break;
           }
           // Correct the ordering within each group ...
           if(n_in_group>1){
             // ... sort halo scores (ascending) ...
             merge_sort(halo_score,      (size_t)n_in_group,&halo_score_index,SID_INT,   SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
             merge_sort(halo_score_index,(size_t)n_in_group,&halo_score_rank, SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
             // ... set the new pointers ...
             new_first=neighbours[group_ids_index[i_halo+halo_score_index[n_in_group-1]]];
             new_last =neighbours[group_ids_index[i_halo+halo_score_index[0]]];
             for(j_halo=i_halo,k_halo=0;k_halo<n_in_group;j_halo++,k_halo++){
               if(neighbours[group_ids_index[j_halo]]!=new_last)
                 neighbours[group_ids_index[j_halo]]->group_halo_next=neighbours[group_ids_index[i_halo+halo_score_index[halo_score_rank[k_halo]-1]]];
               else
                 neighbours[group_ids_index[j_halo]]->group_halo_next=NULL;
               neighbours[group_ids_index[j_halo]]->group_halo_first=new_first;
             }
             SID_free(SID_FARG halo_score_index);
             SID_free(SID_FARG halo_score_rank);
           }
           i_halo+=n_in_group;
         }
         SID_free(SID_FARG halo_score);
       }

       SID_free(SID_FARG group_ids);
       SID_free(SID_FARG group_ids_index);
       SID_free(SID_FARG neighbours);
     }
  }
}
*/

