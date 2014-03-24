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

void read_trees_pointers(tree_info        *trees,
                         const char       *filename_input_dir_horizontal_trees,
                         int               i_file_ptrs,
                         int               i_read_ptrs,
                         tree_node_info ***pointers_groups_local,
                         tree_node_info ***pointers_subgroups_local,
                         int               mode){

   char pointer_type_text[32];
   int  flag_bridge_forematch=FALSE;
   int  flag_bridge_backmatch=FALSE;
   if(check_mode_for_flag(mode,READ_TREES_POINTERS_BRIDGE_FOREMATCH)){
      sprintf(pointer_type_text,"bridge_forematch");
      flag_bridge_forematch=TRUE;
   }
   else if(check_mode_for_flag(mode,READ_TREES_POINTERS_BRIDGE_BACKMATCH)){
      sprintf(pointer_type_text,"bridge_backmatch");
      flag_bridge_backmatch=TRUE;
   }
   else
      SID_trap_error("Invalid mode (%d) passed to read_trees_pointers().",ERROR_LOGIC,mode);

   SID_log("Reading %s pointers for snapshot %03d...",SID_LOG_OPEN,pointer_type_text,trees->snap_list[i_file_ptrs]);

   // Allocate array
   pointers_groups_local[i_file_ptrs]=
      (tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_groups_snap_local[i_file_ptrs]);
   pointers_subgroups_local[i_file_ptrs]=
      (tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_subgroups_snap_local[i_file_ptrs]);

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
   sprintf(filename_ptrs_in,"%s/horizontal_trees_%s_pointers_%03d.dat",filename_input_dir_horizontal_trees,pointer_type_text,i_read_ptrs);
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
     int group_tree_id;
     int group_snap;
     int group_index;
     int n_subgroups_group;
     SID_fread_all(tree_read_buffer,4*sizeof(int),1,&fp_ptrs_in);
     group_tree_id    =tree_read_buffer[0];
     group_snap       =tree_read_buffer[1];
     group_index      =tree_read_buffer[2];
     n_subgroups_group=tree_read_buffer[3];
     
     // Read each subgroup in turn
     for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
       int subgroup_tree_id;
       int subgroup_snap;
       int subgroup_index;
       SID_fread_all(tree_read_buffer,3*sizeof(int),1,&fp_ptrs_in);
       subgroup_tree_id=tree_read_buffer[0];
       subgroup_snap   =tree_read_buffer[1];
       subgroup_index  =tree_read_buffer[2];

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

            // Sanity check
            if(flag_bridge_forematch){
               if(check_mode_for_flag(result_i->tree_case,TREE_CASE_MATCHED_TO_BRIDGE) && result_j==NULL)
                  SID_trap_error("Could not find the bridge forematch a group (%d,%d) marked TREE_CASE_MATCHED_TO_BRIDGE.",ERROR_LOGIC,i_file_ptrs,i_group);
            }
            else if(flag_bridge_backmatch){
               if(check_mode_for_flag(result_i->tree_case,TREE_CASE_EMERGED_CANDIDATE) && result_j==NULL)
                  SID_trap_error("Could not find the bridge forematch a group (%d,%d) marked TREE_CASE_MATCHED_TO_BRIDGE.",ERROR_LOGIC,i_file_ptrs,i_group);
            }

            // Create pointer
            pointers_groups_local[result_i->snap_tree][result_i->neighbour_index]=result_j;

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

         // Sanity check
         if(flag_bridge_forematch){
            if(check_mode_for_flag(result_i->tree_case,TREE_CASE_MATCHED_TO_BRIDGE) && result_j==NULL)
              SID_trap_error("Could not find a subgroup bridge forematch (%d,%d).",ERROR_LOGIC,i_file_ptrs,i_subgroup);
         }
         else if(flag_bridge_backmatch){
            if(check_mode_for_flag(result_i->tree_case,TREE_CASE_EMERGED_CANDIDATE) && result_j==NULL)
              SID_trap_error("Could not find a subgroup bridge backmatch (%d,%d).",ERROR_LOGIC,i_file_ptrs,i_subgroup);
         }

         // Create pointer
         pointers_subgroups_local[result_i->snap_tree][result_i->neighbour_index]=result_j;

         current_subgroup_local=current_subgroup_local->next_neighbour;
       }
     }
   } // i_group
   SID_fclose(&fp_ptrs_in);

   SID_log("Done.",SID_LOG_CLOSE);
}

