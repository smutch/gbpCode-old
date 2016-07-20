#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_fragmented_info(tree_horizontal_extended_info **groups,   int *n_groups,
                               tree_horizontal_extended_info **subgroups,int *n_subgroups,
                               int        **n_subgroups_group,
                               int          i_read, // tree snapshot index
                               int          j_read, // actual snapshot index
                               int          l_read,
                               int          i_read_step,
                               int          n_wrap){
   SID_log("Propagating fragmented halo information for snapshot #%03d...",SID_LOG_OPEN,j_read);
   // Process groups
   int i_group;
   int i_subgroup;
   int j_subgroup;
   int flag_returned;
   for(i_group=0,i_subgroup=0;i_group<n_groups[l_read];i_group++){
      int group_id                =groups[i_read%n_wrap][i_group].id;
      int group_tree_id           =groups[i_read%n_wrap][i_group].tree_id;
      int group_descendant_id     =groups[i_read%n_wrap][i_group].descendant_id;
      int group_type              =groups[i_read%n_wrap][i_group].type;
      int group_file_offset       =groups[i_read%n_wrap][i_group].descendant_file_offset;
      int group_index             =groups[i_read%n_wrap][i_group].descendant_index;
      int group_score_desc        =groups[i_read%n_wrap][i_group].score_desc;
      int group_score_prog        =groups[i_read%n_wrap][i_group].score_prog;
      int group_snap_bridge       =groups[i_read%n_wrap][i_group].snap_bridge;
      int group_file_bridge       =groups[i_read%n_wrap][i_group].file_bridge;
      int group_index_bridge      =groups[i_read%n_wrap][i_group].index_bridge;
      int group_id_bridge         =groups[i_read%n_wrap][i_group].id_bridge;

      // Check if the group's descendant has returned to it's bridge (if fragmented)
      flag_returned=(group_file_bridge==(i_read+group_file_offset) && group_index_bridge==group_index);

      // Propagate type.  Stop when the fragmented halo 
      //    returns to the main progenitor line of it's bridge.
      if(group_id==group_descendant_id && !flag_returned){
         if(group_index>=0){ // Important for strayed cases
            // Add the propagated type.
            if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_STRAYED))
               groups[(i_read+group_file_offset)%n_wrap][group_index].type|=TREE_CASE_FRAGMENTED_STRAYED;
            if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_NORMAL))
               groups[(i_read+group_file_offset)%n_wrap][group_index].type|=TREE_CASE_FRAGMENTED_NORMAL;
            if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_OTHER))
               groups[(i_read+group_file_offset)%n_wrap][group_index].type|=TREE_CASE_FRAGMENTED_OTHER;
            // Count the number of flags the descendant already has switched on
            int i_count=0;
            if(check_mode_for_flag(groups[(i_read+group_file_offset)%n_wrap][group_index].type,TREE_CASE_FRAGMENTED_STRAYED)) i_count++;
            if(check_mode_for_flag(groups[(i_read+group_file_offset)%n_wrap][group_index].type,TREE_CASE_FRAGMENTED_NORMAL))  i_count++;
            if(check_mode_for_flag(groups[(i_read+group_file_offset)%n_wrap][group_index].type,TREE_CASE_FRAGMENTED_OTHER)) i_count++;
            // Check that the affected halo does not have more than one fragmented halo flag turned on
            if(i_count>1)
               SID_trap_error("Multiple (%d) TREE_CASE_FRAGMENT switches present (type=%d) for i_snap/i_group=%d/%d when a max of one is alowed. Progenitor info: i_snap/i_halo/type=%d/%d/%d",
                              ERROR_LOGIC,
                              i_count,
                              groups[(i_read+group_file_offset)%n_wrap][group_index].type,
                              j_read+group_file_offset*i_read_step,
                              group_index,
                              j_read,
                              i_group,
                              group_type);

         }
      }

      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         int subgroup_id                =subgroups[i_read%n_wrap][i_subgroup].id;
         int subgroup_tree_id           =subgroups[i_read%n_wrap][i_subgroup].tree_id;
         int subgroup_descendant_id     =subgroups[i_read%n_wrap][i_subgroup].descendant_id;
         int subgroup_type              =subgroups[i_read%n_wrap][i_subgroup].type;
         int subgroup_score_desc        =subgroups[i_read%n_wrap][i_subgroup].score_desc;
         int subgroup_score_prog        =subgroups[i_read%n_wrap][i_subgroup].score_prog;
         int subgroup_snap_bridge       =subgroups[i_read%n_wrap][i_subgroup].snap_bridge;
         int subgroup_file_bridge       =subgroups[i_read%n_wrap][i_subgroup].file_bridge;
         int subgroup_index_bridge      =subgroups[i_read%n_wrap][i_subgroup].index_bridge;
         int subgroup_id_bridge         =subgroups[i_read%n_wrap][i_subgroup].id_bridge;
         int subgroup_file_offset       =subgroups[i_read%n_wrap][i_subgroup].descendant_file_offset;
         int subgroup_index             =subgroups[i_read%n_wrap][i_subgroup].descendant_index;

         // Check if the subgroup's descendant has returned to it's bridge (if fragmented)
         flag_returned=(subgroup_file_bridge==(i_read+subgroup_file_offset) && subgroup_index_bridge==subgroup_index);

         // Propagate type.  Stop when the fragmented halo  
         //    returns to the main progenitor line of it's bridge.
         if(subgroup_id==subgroup_descendant_id && !flag_returned){
            if(subgroup_index>=0){ // Important for strayed cases
               // Add the propagated type.
               if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_STRAYED))
                  subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type|=TREE_CASE_FRAGMENTED_STRAYED;
               if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_NORMAL))
                  subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type|=TREE_CASE_FRAGMENTED_NORMAL;
               if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_OTHER))
                  subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type|=TREE_CASE_FRAGMENTED_OTHER;
               // Count the number of flags the descendant already has switched on
               int i_count=0;
               if(check_mode_for_flag(subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,TREE_CASE_FRAGMENTED_STRAYED)) i_count++;
               if(check_mode_for_flag(subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,TREE_CASE_FRAGMENTED_NORMAL))  i_count++;
               if(check_mode_for_flag(subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,TREE_CASE_FRAGMENTED_OTHER)) i_count++;
               // Check that the affected halo does not have more than one fragmented halo flag turned on
               if(i_count>1)
                  SID_trap_error("Multiple (%d) TREE_CASE_FRAGMENT switches present (type=%d) for i_snap/i_subgroup=%d/%d when a max of one is alowed. Progenitor info: i_snap/i_halo/type=%d/%d/%d",
                                 ERROR_LOGIC,
                                 i_count,
                                 subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,
                                 j_read+subgroup_file_offset*i_read_step,
                                 subgroup_index,
                                 j_read,
                                 i_subgroup,
                                 subgroup_type);
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

