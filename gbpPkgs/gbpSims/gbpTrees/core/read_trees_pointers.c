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
                         int               mode){

   char pointer_type_text[32];
   int  flag_bridge_forematch=FALSE;
   int  flag_bridge_backmatch=FALSE;
   tree_node_info ***pointers_groups_local     =NULL;
   tree_node_info ***pointers_subgroups_local  =NULL;
   float           **score_groups_local        =NULL;
   float           **score_groups_prog_local   =NULL;
   float           **score_subgroups_local     =NULL;
   float           **score_subgroups_prog_local=NULL;
   if(check_mode_for_flag(mode,READ_TREES_POINTERS_BRIDGE_FOREMATCH) && !check_mode_for_flag(mode,READ_TREES_POINTERS_BRIDGE_BACKMATCH)){
      sprintf(pointer_type_text,"forematch");
      flag_bridge_forematch=TRUE;
      // Fetch the pointers from the trees structure
      if(trees->group_forematch_pointers==NULL)
         SID_trap_error("Group forematch pointer arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         pointers_groups_local=trees->group_forematch_pointers;
      if(trees->subgroup_forematch_pointers==NULL)
         SID_trap_error("Subgroup forematch pointer arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         pointers_subgroups_local=trees->subgroup_forematch_pointers;
      if(trees->group_forematch_score==NULL)
         SID_trap_error("Group forematch score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_groups_local=trees->group_forematch_score;
      if(trees->subgroup_forematch_score==NULL)
         SID_trap_error("Subgroup forematch score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_subgroups_local=trees->subgroup_forematch_score;
      if(trees->group_descendant_score==NULL)
         SID_trap_error("Group descendant score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_groups_prog_local=trees->group_descendant_score;
      if(trees->subgroup_descendant_score==NULL)
         SID_trap_error("Subgroup descendant score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_subgroups_prog_local=trees->subgroup_descendant_score;
   }
   else if(check_mode_for_flag(mode,READ_TREES_POINTERS_BRIDGE_BACKMATCH) && !check_mode_for_flag(mode,READ_TREES_POINTERS_BRIDGE_FOREMATCH)){
      sprintf(pointer_type_text,"backmatch");
      flag_bridge_backmatch=TRUE;
      // Fetch the pointers from the trees structure
      if(trees->group_backmatch_pointers==NULL)
         SID_trap_error("Group backmatch pointer arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         pointers_groups_local=trees->group_backmatch_pointers;
      if(trees->subgroup_backmatch_pointers==NULL)
         SID_trap_error("Subgroup backmatch pointer arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         pointers_subgroups_local=trees->subgroup_backmatch_pointers;
      if(trees->group_backmatch_score==NULL)
         SID_trap_error("Group backmatch score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_groups_local=trees->group_backmatch_score;
      if(trees->subgroup_backmatch_score==NULL)
         SID_trap_error("Subgroup backmatch score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_subgroups_local=trees->subgroup_backmatch_score;
      if(trees->group_progenitor_score==NULL)
         SID_trap_error("Group progenitor score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_groups_prog_local=trees->group_progenitor_score;
      if(trees->subgroup_progenitor_score==NULL)
         SID_trap_error("Subgroup progenitor score arrays have not been initialized properly before a call to read_trees_pointers().",ERROR_LOGIC);
      else
         score_subgroups_prog_local=trees->subgroup_progenitor_score;
   }
   else
      SID_trap_error("Invalid mode (%d) passed to read_trees_pointers().",ERROR_LOGIC,mode);

   SID_log("Reading %s pointers for snapshot %03d...",SID_LOG_OPEN,pointer_type_text,trees->snap_list[i_file_ptrs]);

   // Allocate array
   pointers_groups_local[i_file_ptrs]=
      (tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_groups_snap_local[i_file_ptrs]);
   pointers_subgroups_local[i_file_ptrs]=
      (tree_node_info **)SID_malloc(sizeof(tree_node_info *)*trees->n_subgroups_snap_local[i_file_ptrs]);
   score_groups_local[i_file_ptrs]=
      (float *)SID_malloc(sizeof(float)*trees->n_groups_snap_local[i_file_ptrs]);
   score_groups_prog_local[i_file_ptrs]=
      (float *)SID_malloc(sizeof(float)*trees->n_groups_snap_local[i_file_ptrs]);
   score_subgroups_local[i_file_ptrs]=
      (float *)SID_malloc(sizeof(float)*trees->n_subgroups_snap_local[i_file_ptrs]);
   score_subgroups_prog_local[i_file_ptrs]=
      (float *)SID_malloc(sizeof(float)*trees->n_subgroups_snap_local[i_file_ptrs]);

   // Open file and read header
   char   filename_ptrs_in[MAX_FILENAME_LENGTH];
   int    n_buffer_group   =4*sizeof(int)+2*sizeof(float);
   int    n_buffer_subgroup=3*sizeof(int)+2*sizeof(float);
   int    n_buffer_alloc   =MAX(n_buffer_group,n_buffer_subgroup);
   char  *tree_read_buffer =(char *)SID_malloc(n_buffer_alloc);
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
   int i_buffer;
   tree_node_info *current_group_local   =trees->first_neighbour_groups[i_file_ptrs];
   tree_node_info *current_subgroup_local=trees->first_neighbour_subgroups[i_file_ptrs];
   for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
     int flag_group_added=FALSE;
     // Read horizontal trees for groups
     int   group_tree_id;
     int   group_snap;
     int   group_index;
     float group_score;
     float group_score_prog;
     int   n_subgroups_group;
     SID_fread_all(tree_read_buffer,n_buffer_group,1,&fp_ptrs_in);
     i_buffer=0;
     group_tree_id    =((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
     group_snap       =((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
     group_index      =((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
     group_score      =((float *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(float);
     group_score_prog =((float *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(float);
     n_subgroups_group=((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
     
     // Read each subgroup in turn
     for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
       int   subgroup_tree_id;
       int   subgroup_snap;
       int   subgroup_index;
       float subgroup_score;
       float subgroup_score_prog;
       SID_fread_all(tree_read_buffer,n_buffer_subgroup,1,&fp_ptrs_in);
       i_buffer=0;
       subgroup_tree_id   =((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
       subgroup_snap      =((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
       subgroup_index     =((int   *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(int);
       subgroup_score     =((float *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(float);
       subgroup_score_prog=((float *)(&(tree_read_buffer[i_buffer])))[0];i_buffer+=sizeof(float);

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

            // Create pointer and store scpre
            pointers_groups_local[result_i->snap_tree][result_i->neighbour_index]=result_j;
            score_groups_local[result_i->snap_tree][result_i->neighbour_index]     =group_score;
            score_groups_prog_local[result_i->snap_tree][result_i->neighbour_index]=group_score_prog;

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

         // Create pointer and store score
         pointers_subgroups_local[result_i->snap_tree][result_i->neighbour_index]=result_j;
         score_subgroups_local[result_i->snap_tree][result_i->neighbour_index]     =subgroup_score;
         score_subgroups_prog_local[result_i->snap_tree][result_i->neighbour_index]=subgroup_score_prog;

         current_subgroup_local=current_subgroup_local->next_neighbour;
       }
     }
   } // i_group
   SID_free(SID_FARG tree_read_buffer);
   SID_fclose(&fp_ptrs_in);

   SID_log("Done.",SID_LOG_CLOSE);
}

