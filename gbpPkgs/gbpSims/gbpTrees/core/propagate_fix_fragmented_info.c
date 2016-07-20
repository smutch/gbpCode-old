#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_fix_fragmented_info(tree_horizontal_extended_info **groups,   int *n_groups,
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

      // For every halo that is fragmented in this snapshot, we need
      //    to check if it should be.  Scan their progenitors for a merger primary.  
      //    If it isn't fragmented, then this halo shouldn't be either.  If it is, 
      //    turn-off fragmented flags and recompute the peak particle count.
      if(check_if_type_is_fragmented(group_type)){
         tree_horizontal_extended_info *group        =&(groups[i_read%n_wrap][i_group]);
         tree_horizontal_extended_info *current_group=set_extended_first_progenitor(groups,group,n_wrap);
         int flag_fix=FALSE;
         int n_p_peak=0; // calculated below under the assumption that a fix is being made
         while(current_group!=NULL){
            int current_group_type=current_group->type;
            if(!check_if_type_is_fragmented(current_group_type)){
               n_p_peak=MAX(n_p_peak,current_group->n_particles_peak);
               if(check_mode_for_flag(current_group_type,TREE_CASE_MERGER_PRIMARY))
                  flag_fix=TRUE;
            }
            current_group=set_extended_next_progenitor(groups,current_group,n_wrap);
         }
         if(flag_fix){
            // Apply dominant halo condition to complete calculation of peak particle count
            if(check_mode_for_flag(group_type,TREE_CASE_DOMINANT)==check_mode_for_flag(group_type,TREE_CASE_MOST_MASSIVE))
               n_p_peak=MAX(n_p_peak,group->n_particles);
            // Apply fix
            group->n_particles_peak =n_p_peak;
            group->type            &=(~TREE_CASE_FRAGMENTED_STRAYED);
            group->type            &=(~TREE_CASE_FRAGMENTED_NORMAL);
            group->type            &=(~TREE_CASE_FRAGMENTED_OTHER);
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

         // For every halo that is fragmented in this snapshot, we need
         //    to check if it should be.  Scan their progenitors for a merger primary.
         //    If it isn't fragmented, then this halo shouldn't be either.  If it is,
         //    turn-off fragmented flags and recompute the peak particle count.
         if(check_if_type_is_fragmented(subgroup_type)){
            tree_horizontal_extended_info *subgroup        =&(subgroups[i_read%n_wrap][i_subgroup]);
            tree_horizontal_extended_info *current_subgroup=set_extended_first_progenitor(subgroups,subgroup,n_wrap);
            int flag_fix=FALSE;
            int n_p_peak=0; // calculated below under the assumption that a fix is being made
            while(current_subgroup!=NULL){
               int current_subgroup_type=current_subgroup->type;
               if(!check_if_type_is_fragmented(current_subgroup_type)){
                  n_p_peak=MAX(n_p_peak,current_subgroup->n_particles_peak);
                  if(check_mode_for_flag(current_subgroup_type,TREE_CASE_MERGER_PRIMARY))
                     flag_fix=TRUE;
               }
               current_subgroup=set_extended_next_progenitor(subgroups,current_subgroup,n_wrap);
            }
            if(flag_fix){
               // Apply dominant halo condition to complete calculation of peak particle count
               if(check_mode_for_flag(subgroup_type,TREE_CASE_DOMINANT)==check_mode_for_flag(subgroup_type,TREE_CASE_MOST_MASSIVE))
                  n_p_peak=MAX(n_p_peak,subgroup->n_particles);
               // Apply fix
               subgroup->n_particles_peak =n_p_peak;
               subgroup->type            &=(~TREE_CASE_FRAGMENTED_STRAYED);
               subgroup->type            &=(~TREE_CASE_FRAGMENTED_NORMAL);
               subgroup->type            &=(~TREE_CASE_FRAGMENTED_OTHER);
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

