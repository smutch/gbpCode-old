#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpZFIRE.h>

void analyze_halos_and_N_subhalos(tree_info  *trees,
                                  const char *filename_out_root,
                                  const char *catalog_root,
                                  double      z_obs_exact,
                                  double      M_cut_lo,
                                  double      M_cut_hi,
                                  int         n_subgroups_track_max){

  // Compute merger rates ...
  SID_log("Constructing catalogs...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Fetch halo properties
  halo_properties_info **group_properties   =(halo_properties_info **)ADaPS_fetch(trees->data,"properties_groups");
  halo_properties_info **subgroup_properties=(halo_properties_info **)ADaPS_fetch(trees->data,"properties_subgroups");

  // Determine the best z_obs snapshot
  int    i_z_obs=find_treesnap_z(trees,z_obs_exact);
  int    i_z_0  =trees->n_snaps-1;
  double z_obs  =trees->z_list[i_z_obs];
  double t_obs  =trees->t_list[i_z_obs];
  SID_log("The snapshot corresponding best to z=%4.2f is #%03d (z=%4.2f).",SID_LOG_COMMENT,z_obs_exact,trees->snap_list[i_z_obs],z_obs);

  // Allocate the list arrays
  tree_node_info  **list_groups       =(tree_node_info  **)SID_malloc(sizeof(tree_node_info  *)*trees->n_groups_snap_local[i_z_obs]);
  tree_node_info  **list_subgroups_all=(tree_node_info  **)SID_malloc(sizeof(tree_node_info  *)*trees->n_groups_snap_local[i_z_obs]);
  tree_node_info ***list_subgroups    =(tree_node_info ***)SID_malloc(sizeof(tree_node_info **)*n_subgroups_track_max);
  for(int i_track=0;i_track<n_subgroups_track_max;i_track++)
     list_subgroups[i_track]=(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_groups_snap_local[i_z_obs]);

  // Select z_obs systems
  tree_node_info *current_group;
  int  n_list_groups       =0;
  int  n_list_subgroups_all=0;
  int *n_list_subgroups    =(int *)SID_calloc(sizeof(int)*n_subgroups_track_max);
  SID_log("Selecting z=%4.2f systems...",SID_LOG_OPEN,trees->z_list[i_z_obs]);
  current_group=trees->first_neighbour_groups[i_z_obs];
  while(current_group!=NULL){
     halo_properties_info *current_group_properties=&(group_properties[current_group->snap_tree][current_group->neighbour_index]);
     if(current_group_properties->M_vir>=M_cut_lo && current_group_properties->M_vir<=M_cut_hi){
        list_groups[n_list_groups++]=current_group;
        tree_node_info *current_subgroup=current_group->substructure_first;
        int i_subgroup       =0;
        int i_subgroups_track=0;
        while(current_subgroup!=NULL && i_subgroups_track<n_subgroups_track_max){
           if(!check_treenode_if_fragmented(current_subgroup)){
              if(i_subgroup>0)
                 list_subgroups_all[n_list_subgroups_all++]=current_subgroup;
              list_subgroups[i_subgroups_track][n_list_subgroups[i_subgroups_track]++]=current_subgroup;
              i_subgroups_track++;
           }
           i_subgroup++;
           current_subgroup=current_subgroup->substructure_next;
        }
     }
     current_group=current_group->next_neighbour;
  }
  SID_log("%d groups found...",SID_LOG_CONTINUE,n_list_groups);
  SID_log("Done.",SID_LOG_CLOSE);

  // Write properties and tracks for selected z_obs groups and subgroups
  char catalog_name[32];
  sprintf(catalog_name,"%s_groups",catalog_root);
  write_tree_branches(trees,list_groups,n_list_groups,TRUE,filename_out_root,catalog_name);
  average_tree_branches(catalog_name);
  for(int i_track=0;i_track<n_subgroups_track_max;i_track++){
     if(i_track==0){
        sprintf(catalog_name,"%s_subgroups_all",catalog_root);
        write_tree_branches(trees,list_subgroups_all,n_list_subgroups_all,FALSE,filename_out_root,catalog_name);
        average_tree_branches(catalog_name);
     }
     sprintf(catalog_name,"%s_subgroups_%02d",catalog_root,i_track);
     write_tree_branches(trees,list_subgroups[i_track],n_list_subgroups[i_track],FALSE,filename_out_root,catalog_name);
     average_tree_branches(catalog_name);
  }

  // Clean-up
  SID_free(SID_FARG n_list_subgroups);
  for(int i_track=0;i_track<n_subgroups_track_max;i_track++)
     SID_free(SID_FARG list_subgroups[i_track]);
  SID_free(SID_FARG list_subgroups);
  SID_free(SID_FARG list_subgroups_all);
  SID_free(SID_FARG list_groups);

  SID_log("Done.",SID_LOG_CLOSE);
}

