#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_bridge_info(tree_horizontal_extended_info **groups,   int *n_groups,
                           tree_horizontal_extended_info **subgroups,int *n_subgroups,
                           int        **n_subgroups_group,
                           int          i_read, // tree snapshot index
                           int          j_read, // actual snapshot index
                           int          l_read,
                           int          i_read_step,
                           int          n_wrap){
   SID_log("Propagating bridged halo information for snapshot #%03d...",SID_LOG_OPEN,j_read);
   // Process groups
   int i_group;
   int i_subgroup;
   int j_subgroup;
   int flag_returned;
   for(i_group=0,i_subgroup=0;i_group<n_groups[l_read];i_group++){
      int group_type              =groups[i_read%n_wrap][i_group].type;
      int group_file_offset       =groups[i_read%n_wrap][i_group].descendant_file_offset;
      int group_index             =groups[i_read%n_wrap][i_group].descendant_index;
      int group_snap_bridge       =groups[i_read%n_wrap][i_group].snap_bridge;
      int group_file_bridge       =groups[i_read%n_wrap][i_group].file_bridge;
      int group_index_bridge      =groups[i_read%n_wrap][i_group].index_bridge;
      int group_id_bridge         =groups[i_read%n_wrap][i_group].id_bridge;

      // Advance the bridge to it's current-time descendant. Note that n_wrap needs to be
      //   at least 2*n_search+2 so that we are sure to have the bridged halo as well
      while(group_file_bridge<=i_read && group_index_bridge>=0){
         tree_horizontal_extended_info *bridge_descendant=&(groups[group_file_bridge%n_wrap][group_index_bridge]);
         group_snap_bridge+=bridge_descendant->descendant_file_offset*i_read_step;
         group_file_bridge+=bridge_descendant->descendant_file_offset;
         group_index_bridge=bridge_descendant->descendant_index;
      }

      // Propagate the information 
      if((check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_STRAYED)  ||
          check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_RETURNED) ||
          check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_EXCHANGED)) &&
          check_mode_for_flag(group_type,TREE_CASE_MAIN_PROGENITOR)){
         if(group_index>=0){ // Important for strayed cases
            int desc_file;
            desc_file=i_read+group_file_offset;
            groups[desc_file%n_wrap][group_index].id_bridge   =group_id_bridge;
            groups[desc_file%n_wrap][group_index].snap_bridge =group_snap_bridge;
            groups[desc_file%n_wrap][group_index].file_bridge =group_file_bridge;
            groups[desc_file%n_wrap][group_index].index_bridge=group_index_bridge;
         }
      }

      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         int subgroup_type              =subgroups[i_read%n_wrap][i_subgroup].type;
         int subgroup_file_offset       =subgroups[i_read%n_wrap][i_subgroup].descendant_file_offset;
         int subgroup_index             =subgroups[i_read%n_wrap][i_subgroup].descendant_index;
         int subgroup_snap_bridge       =subgroups[i_read%n_wrap][i_subgroup].snap_bridge;
         int subgroup_file_bridge       =subgroups[i_read%n_wrap][i_subgroup].file_bridge;
         int subgroup_index_bridge      =subgroups[i_read%n_wrap][i_subgroup].index_bridge;
         int subgroup_id_bridge         =subgroups[i_read%n_wrap][i_subgroup].id_bridge;

         // Advance the bridge to it's current-time descendant. Note that n_wrap needs to be
         //   at least 2*n_search+2 so that we are sure to have the bridged halo as well
         while(subgroup_file_bridge<i_read && subgroup_index_bridge>=0){
            tree_horizontal_extended_info *bridge_descendant=&(subgroups[subgroup_file_bridge%n_wrap][subgroup_index_bridge]);
            subgroup_snap_bridge+=bridge_descendant->descendant_file_offset*i_read_step;
            subgroup_file_bridge+=bridge_descendant->descendant_file_offset;
            subgroup_index_bridge=bridge_descendant->descendant_index;
         }

         // Propagate the information
         if((check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_STRAYED)  ||
             check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_RETURNED) ||
             check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_EXCHANGED)) &&
             check_mode_for_flag(subgroup_type,TREE_CASE_MAIN_PROGENITOR)){
            if(subgroup_index>=0){ // Important for strayed cases
               int desc_file;
               desc_file=i_read+group_file_offset;
               subgroups[desc_file%n_wrap][subgroup_index].id_bridge   =subgroup_id_bridge;
               subgroups[desc_file%n_wrap][subgroup_index].snap_bridge =subgroup_snap_bridge;
               subgroups[desc_file%n_wrap][subgroup_index].file_bridge =subgroup_file_bridge;
               subgroups[desc_file%n_wrap][subgroup_index].index_bridge=subgroup_index_bridge;
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

