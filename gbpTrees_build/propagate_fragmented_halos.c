#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_fragmented_halos(tree_horizontal_extended_info **groups,   int *n_groups,
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
   for(i_group=0,i_subgroup=0;i_group<n_groups[l_read];i_group++){
      int group_id;
      int group_n_particles;
      int group_tree_id;
      int group_descendant_id;
      int group_type;
      int group_type_descendant;
      int group_file_offset;
      int group_n_particles_parent;
      int group_n_particles_desc;
      int group_n_particles_proj;
      int group_score_desc;
      int group_score_prog;
      int group_snap_bridge;
      int group_index_bridge;
      int group_id_bridge;
      int group_index;
      group_id                =groups[i_read%n_wrap][i_group].id;
      group_n_particles       =groups[i_read%n_wrap][i_group].n_particles;
      group_tree_id           =groups[i_read%n_wrap][i_group].tree_id;
      group_descendant_id     =groups[i_read%n_wrap][i_group].descendant_id;
      group_type              =groups[i_read%n_wrap][i_group].type;
      group_file_offset       =groups[i_read%n_wrap][i_group].file_offset;
      group_n_particles_parent=groups[i_read%n_wrap][i_group].n_particles_parent;
      group_n_particles_desc  =groups[i_read%n_wrap][i_group].n_particles_desc;
      group_n_particles_proj  =groups[i_read%n_wrap][i_group].n_particles_proj;
      group_score_desc        =groups[i_read%n_wrap][i_group].score_desc;
      group_score_prog        =groups[i_read%n_wrap][i_group].score_prog;
      group_snap_bridge       =groups[i_read%n_wrap][i_group].snap_bridge;
      group_index_bridge      =groups[i_read%n_wrap][i_group].index_bridge;
      group_id_bridge         =groups[i_read%n_wrap][i_group].id_bridge;
      group_index             =groups[i_read%n_wrap][i_group].index;

      // Propagate type.  This check is needed so that we don't
      //   propagate this type into any halos this one merges with.
      if(group_id==group_descendant_id){
         if(group_index>=0){ // Important for strayed cases
            // Add the propagated type.  The check against TREE_CASE_MERGER is also
            //    needed so that we don't propagate this type into any halos this one merges with.
            if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_STRAYED)   && !check_mode_for_flag(group_type,TREE_CASE_MERGER))
               groups[(i_read+group_file_offset)%n_wrap][group_index].type|=TREE_CASE_FRAGMENTED_STRAYED;
            if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_RETURNED)  && !check_mode_for_flag(group_type,TREE_CASE_MERGER))
               groups[(i_read+group_file_offset)%n_wrap][group_index].type|=TREE_CASE_FRAGMENTED_RETURNED;
            if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_EXCHANGED) && !check_mode_for_flag(group_type,TREE_CASE_MERGER))
               groups[(i_read+group_file_offset)%n_wrap][group_index].type|=TREE_CASE_FRAGMENTED_EXCHANGED;
            // Count the number of flags the descendant already has switched on
            int i_count=0;
            if(check_mode_for_flag(groups[(i_read+group_file_offset)%n_wrap][group_index].type,TREE_CASE_FRAGMENTED_STRAYED))   i_count++;
            if(check_mode_for_flag(groups[(i_read+group_file_offset)%n_wrap][group_index].type,TREE_CASE_FRAGMENTED_RETURNED))  i_count++;
            if(check_mode_for_flag(groups[(i_read+group_file_offset)%n_wrap][group_index].type,TREE_CASE_FRAGMENTED_EXCHANGED)) i_count++;
            // Check that the affected halo does not have more than one fragmented halo flags turned on
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
      /*
      // Don't count fragmented halos as mergers when they're accreated
      else if(check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_STRAYED)  ||
              check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_RETURNED) ||
              check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_EXCHANGED)){
         groups[i_read%n_wrap][i_group].type&=(~TREE_CASE_MERGER);
      }
      */

      // Propagate the information about the bridge as well
      if((check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_STRAYED)  ||
          check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_RETURNED) ||
          check_mode_for_flag(group_type,TREE_CASE_FRAGMENTED_EXCHANGED)) &&
         !check_mode_for_flag(group_type,TREE_CASE_MERGER)){
         if(group_index>=0){ // Important for strayed cases
            int desc_file;
            desc_file=i_read+group_file_offset;
            groups[desc_file%n_wrap][group_index].id_bridge   =group_id_bridge;
            groups[desc_file%n_wrap][group_index].snap_bridge =group_snap_bridge;
            groups[desc_file%n_wrap][group_index].index_bridge=group_index_bridge;
         }
      }

      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         int subgroup_id;
         int subgroup_n_particles;
         int subgroup_tree_id;
         int subgroup_descendant_id;
         int subgroup_type;
         int subgroup_file_offset;
         int subgroup_n_particles_parent;
         int subgroup_n_particles_desc;
         int subgroup_n_particles_proj;
         int subgroup_score_desc;
         int subgroup_score_prog;
         int subgroup_snap_bridge;
         int subgroup_index_bridge;
         int subgroup_id_bridge;
         int subgroup_index;
         int subgroup_type_descendant;;
         subgroup_id                =subgroups[i_read%n_wrap][i_subgroup].id;
         subgroup_n_particles       =subgroups[i_read%n_wrap][i_subgroup].n_particles;
         subgroup_tree_id           =subgroups[i_read%n_wrap][i_subgroup].tree_id;
         subgroup_descendant_id     =subgroups[i_read%n_wrap][i_subgroup].descendant_id;
         subgroup_type              =subgroups[i_read%n_wrap][i_subgroup].type;
         subgroup_file_offset       =subgroups[i_read%n_wrap][i_subgroup].file_offset;
         subgroup_n_particles_parent=subgroups[i_read%n_wrap][i_subgroup].n_particles_parent;
         subgroup_n_particles_desc  =subgroups[i_read%n_wrap][i_subgroup].n_particles_desc;
         subgroup_n_particles_proj  =subgroups[i_read%n_wrap][i_subgroup].n_particles_proj;
         subgroup_score_desc        =subgroups[i_read%n_wrap][i_subgroup].score_desc;
         subgroup_score_prog        =subgroups[i_read%n_wrap][i_subgroup].score_prog;
         subgroup_snap_bridge       =subgroups[i_read%n_wrap][i_subgroup].snap_bridge;
         subgroup_index_bridge      =subgroups[i_read%n_wrap][i_subgroup].index_bridge;
         subgroup_id_bridge         =subgroups[i_read%n_wrap][i_subgroup].id_bridge;
         subgroup_index             =subgroups[i_read%n_wrap][i_subgroup].index;

         // Propagate type.  This check is needed so that we don't
         //   propagate this type into any halos this one merges with.
         if(subgroup_id==subgroup_descendant_id){
            if(subgroup_index>=0){ // Important for strayed cases
               // Add the propagated type.  The check against TREE_CASE_MERGER is also
               //    needed so that we don't propagate this type into any halos this one merges with.
               if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_STRAYED)   && !check_mode_for_flag(subgroup_type,TREE_CASE_MERGER))
                  subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type|=TREE_CASE_FRAGMENTED_STRAYED;
               if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_RETURNED)  && !check_mode_for_flag(subgroup_type,TREE_CASE_MERGER))
                  subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type|=TREE_CASE_FRAGMENTED_RETURNED;
               if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_EXCHANGED) && !check_mode_for_flag(subgroup_type,TREE_CASE_MERGER))
                  subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type|=TREE_CASE_FRAGMENTED_EXCHANGED;
               // Count the number of flags the descendant already has switched on
               int i_count=0;
               if(check_mode_for_flag(subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,TREE_CASE_FRAGMENTED_STRAYED))   i_count++;
               if(check_mode_for_flag(subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,TREE_CASE_FRAGMENTED_RETURNED))  i_count++;
               if(check_mode_for_flag(subgroups[(i_read+subgroup_file_offset)%n_wrap][subgroup_index].type,TREE_CASE_FRAGMENTED_EXCHANGED)) i_count++;
               // Check that the affected halo does not have more than one fragmented halo flags turned on
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

         /*
         // Don't count fragmented halos as mergers when they're accreated
         else if(check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_STRAYED)  ||
                 check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_RETURNED) ||
                 check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_EXCHANGED)){
            subgroups[i_read%n_wrap][i_subgroup].type&=(~TREE_CASE_MERGER);
         }
         */

         // Propagate the information about the bridge as well
         if((check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_STRAYED)  ||
             check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_RETURNED) ||
             check_mode_for_flag(subgroup_type,TREE_CASE_FRAGMENTED_EXCHANGED)) &&
            !check_mode_for_flag(subgroup_type,TREE_CASE_MERGER)){
            if(subgroup_index>=0){ // Important for strayed cases
               int desc_file;
               desc_file=i_read+group_file_offset;
               subgroups[desc_file%n_wrap][subgroup_index].id_bridge   =subgroup_id_bridge;
               subgroups[desc_file%n_wrap][subgroup_index].snap_bridge =subgroup_snap_bridge;
               subgroups[desc_file%n_wrap][subgroup_index].index_bridge=subgroup_index_bridge;
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

