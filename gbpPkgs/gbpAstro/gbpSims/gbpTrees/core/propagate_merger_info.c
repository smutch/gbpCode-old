#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_merger_info(tree_horizontal_extended_info **groups,   int *n_groups,
                           tree_horizontal_extended_info **subgroups,int *n_subgroups,
                           int        **n_subgroups_group,
                           int          i_read, // tree snapshot index
                           int          j_read, // actual snapshot index
                           int          l_read,
                           int          i_read_step,
                           int          n_wrap){
   SID_log("Propagating merger information for snapshot #%03d...",SID_LOG_OPEN,j_read);
   // Process groups
   int i_group;
   int i_subgroup;
   int j_subgroup;
   int flag_returned;
   for(i_group=0,i_subgroup=0;i_group<n_groups[l_read];i_group++){
      tree_horizontal_extended_info *this_group=&(groups[i_read%n_wrap][i_group]);
      int group_id                  =this_group->id;
      int group_n_particles         =this_group->n_particles;
      int group_tree_id             =this_group->tree_id;
      int group_descendant_id       =this_group->descendant_id;
      int group_file_offset         =this_group->file_offset;
      int group_index               =this_group->index;
      int group_n_particles_parent  =this_group->n_particles_parent;
      int group_n_particles_desc    =this_group->n_particles_desc;
      int group_n_particles_proj    =this_group->n_particles_proj;
      int group_score_desc          =this_group->score_desc;
      int group_score_prog          =this_group->score_prog;
      int group_snap_bridge         =this_group->snap_bridge;
      int group_file_bridge         =this_group->file_bridge;
      int group_index_bridge        =this_group->index_bridge;
      int group_id_bridge           =this_group->id_bridge;
      tree_horizontal_extended_info *this_desc =&(groups[(i_read+group_file_offset)%n_wrap][group_index]);

      // Set dominant halo flag
/*
      if(check_mode_for_flag(this_group->type,TREE_CASE_NO_PROGENITOR))
         this_group->type|=TREE_CASE_DOMINANT;

      // Propagate dominant halo flags
      if(check_mode_for_flag(this_group->type,TREE_CASE_DOMINANT))
         this_desc->type|=TREE_CASE_DOMINANT;

      // Set merger flags...
      //    ... first, find the primary halo ...
      tree_horizontal_extended_info *primary_group         =this_group;
      int                            primary_group_dominant=check_mode_for_flag(primary_group->type,TREE_CASE_DOMINANT);
      match_info *current_progenitor=&(this_desc->first_progenitor);
      while(current_progenitor->halo!=NULL){
         if(current_progenitor->halo!=this_group){
            int current_progenitor_dominant=check_mode_for_flag(current_progenitor->halo->type,TREE_CASE_DOMINANT);
            int set_new_primary            =FALSE;
            // Decide if this halo is a better primary
            if(current_progenitor_dominant){
               if(!primary_group_dominant)
                  set_new_primary=TRUE;
               else if()
                  set_new_primary=TRUE;
            }
            else if(!primary_group_dominant){
            }
            // ... set the new primary if so.
            if(set_new_primary){
            }
         }
         current_progenitor=current_progenitor->halo.next_progenitor;
      }

      //    ... if this halo is not the primary, set it's merger flag.
      if(this_group!=primary_group)
         this_group->type|=TREE_CASE_MERGER;
*/

      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         int subgroup_id                =subgroups[i_read%n_wrap][i_subgroup].id;
         int subgroup_n_particles       =subgroups[i_read%n_wrap][i_subgroup].n_particles;
         int subgroup_tree_id           =subgroups[i_read%n_wrap][i_subgroup].tree_id;
         int subgroup_descendant_id     =subgroups[i_read%n_wrap][i_subgroup].descendant_id;
         int subgroup_type              =subgroups[i_read%n_wrap][i_subgroup].type;
         int subgroup_file_offset       =subgroups[i_read%n_wrap][i_subgroup].file_offset;
         int subgroup_index             =subgroups[i_read%n_wrap][i_subgroup].index;
         int subgroup_n_particles_parent=subgroups[i_read%n_wrap][i_subgroup].n_particles_parent;
         int subgroup_n_particles_desc  =subgroups[i_read%n_wrap][i_subgroup].n_particles_desc;
         int subgroup_n_particles_proj  =subgroups[i_read%n_wrap][i_subgroup].n_particles_proj;
         int subgroup_score_desc        =subgroups[i_read%n_wrap][i_subgroup].score_desc;
         int subgroup_score_prog        =subgroups[i_read%n_wrap][i_subgroup].score_prog;
         int subgroup_snap_bridge       =subgroups[i_read%n_wrap][i_subgroup].snap_bridge;
         int subgroup_file_bridge       =subgroups[i_read%n_wrap][i_subgroup].file_bridge;
         int subgroup_index_bridge      =subgroups[i_read%n_wrap][i_subgroup].index_bridge;
         int subgroup_id_bridge         =subgroups[i_read%n_wrap][i_subgroup].id_bridge;

      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

