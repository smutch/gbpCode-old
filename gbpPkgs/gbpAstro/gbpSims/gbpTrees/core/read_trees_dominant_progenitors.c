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

void read_trees_dominant_progenitors(tree_info        *trees,
                                     const char       *filename_input_dir_horizontal_trees,
                                     int               i_file_ptrs,
                                     int               i_read_ptrs){

   SID_log("Reading dominanat progenitor pointers for snapshot %03d...",SID_LOG_OPEN,trees->snap_list[i_file_ptrs]);

   // Open file and read header
   char   filename_ptrs_in[MAX_FILENAME_LENGTH];
   int    tree_read_buffer[4];
   SID_fp fp_ptrs_in;
   int    n_step_in;
   int    n_search_in;
   int    n_groups;
   int    n_subgroups;
   int    n_groups_max_in;
   int    n_subgroups_max_in;
   int    n_trees_subgroup;
   int    n_trees_group;
   sprintf(filename_ptrs_in,"%s/horizontal_trees_dominant_progenitor_pointers_%03d.dat",filename_input_dir_horizontal_trees,i_read_ptrs);
   SID_fopen(filename_ptrs_in,"r",&fp_ptrs_in);
   SID_fread_all(&n_step_in,         sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_search_in,       sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_groups,          sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_ptrs_in);
   SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_ptrs_in);

   int i_group;
   int i_subgroup;
   int j_subgroup;
   tree_node_info *current_group_local   =trees->first_neighbour_groups[i_file_ptrs];
   tree_node_info *current_subgroup_local=trees->first_neighbour_subgroups[i_file_ptrs];
   for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
     int flag_group_added=FALSE;
     // Read horizontal trees for groups
     SID_fread_all(tree_read_buffer,4*sizeof(int),1,&fp_ptrs_in);
     int group_tree_id    =tree_read_buffer[0];
     int group_snap       =tree_read_buffer[1];
     int group_index      =tree_read_buffer[2];
     int n_subgroups_group=tree_read_buffer[3];
     
     // Read each subgroup in turn
if(i_group<10) fprintf(stderr,"%d\n",n_subgroups_group);
     for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
       SID_fread_all(tree_read_buffer,3*sizeof(int),1,&fp_ptrs_in);
       int subgroup_tree_id=tree_read_buffer[0];
       int subgroup_snap   =tree_read_buffer[1];
       int subgroup_index  =tree_read_buffer[2];

       // Ignore halos with undefined tree_ids
       int subgroup_forest_id;
       int i_forest;
       if(subgroup_tree_id>=0){
          subgroup_forest_id=trees->tree2forest_mapping_subgroup[subgroup_tree_id];
          i_forest          =subgroup_forest_id-trees->forest_lo_subgroup_local;
       }
       else{
          subgroup_forest_id=-1;
          i_forest          =-1;
       }

       // If this subgroup belongs to a local forest ...
       if(i_forest>=0 && i_forest<trees->n_forests_local){ 
         // ... add the group ...
         tree_node_info *result_i;
         tree_node_info *result_j;
         int             flag_found_j=FALSE;
         if(!flag_group_added){
            // Sanity check
            if(current_group_local->file_index!=i_group)
               SID_trap_error("Group pointer catalog is out of sync (ie. %d!=%d)",ERROR_LOGIC,current_group_local->file_index,i_group);

            // Find the halo pointed to and the halo pointed from (which-is-which depends on the mode)
            result_i    =current_group_local;
            flag_found_j=find_tree_node(trees,group_snap,group_index,TRUE,&result_j);

            // Create pointer
            result_i->progenitor_dominant=result_j;

            flag_group_added=TRUE;
            current_group_local=current_group_local->next_neighbour;
         }

         // ... add the subgroup ...

         // Sanity check
         if(current_subgroup_local->file_index!=i_subgroup)
            SID_trap_error("Subgroup pointer catalog is out of sync (ie. %d!=%d)",ERROR_LOGIC,current_subgroup_local->file_index,i_subgroup);

         // Find the halo pointed to and the halo pointed from (which-is-which depends on the mode)
         result_i    =current_subgroup_local;
         flag_found_j=find_tree_node(trees,subgroup_snap,subgroup_index,FALSE,&result_j);

         // Create pointer
         result_i->progenitor_dominant=result_j;

         current_subgroup_local=current_subgroup_local->next_neighbour;
       }
     }
   } // i_group
   SID_fclose(&fp_ptrs_in);

   SID_log("Done.",SID_LOG_CLOSE);
}

