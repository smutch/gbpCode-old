#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_n_particles_peak(tree_horizontal_extended_info **groups,   int *n_groups,
                                tree_horizontal_extended_info **subgroups,int *n_subgroups,
                                int        **n_subgroups_group,
                                int          i_read, // tree snapshot index
                                int          j_read, // actual snapshot index
                                int          l_read,
                                int          i_read_step,
                                int          n_wrap){
   SID_log("Propagating peak halo sizes for snapshot #%03d...",SID_LOG_OPEN,j_read);
   // Process groups
   int i_subgroup=0;
   for(int i_group=0;i_group<n_groups[l_read];i_group++){
      tree_horizontal_extended_info *this_group=&(groups[i_read%n_wrap][i_group]);

      // Initialize result at leaves
      if(check_mode_for_flag(this_group->type,TREE_CASE_NO_PROGENITORS) || this_group->n_particles_peak==0)
         this_group->n_particles_peak=this_group->n_particles;
      // ... else check current size against a (previously) propagated result.
      else 
         this_group->n_particles_peak=MAX(this_group->n_particles_peak,this_group->n_particles);

      // Propagate peak particle counts forward (unless this halo is a merging fragment)
      tree_horizontal_extended_info *this_group_desc=set_extended_descendant(groups,this_group,i_read,n_wrap);
      if(this_group_desc!=NULL){
         if(check_if_type_is_fragmented(this_group->type)==check_if_type_is_fragmented(this_group_desc->type))
            this_group_desc->n_particles_peak=MAX(this_group_desc->n_particles_peak,this_group->n_particles_peak);
      }

      // Process subgroups
      for(int j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         tree_horizontal_extended_info *this_subgroup=&(subgroups[i_read%n_wrap][i_subgroup]);

         // Initialize peak particle counts at leaves
         int flag_most_massive=check_mode_for_flag(this_subgroup->type,TREE_CASE_MOST_MASSIVE);
         int flag_dominant    =check_mode_for_flag(this_subgroup->type,TREE_CASE_DOMINANT);
         if(check_mode_for_flag(this_subgroup->type,TREE_CASE_NO_PROGENITORS) || this_subgroup->n_particles_peak==0)
            this_subgroup->n_particles_peak=this_subgroup->n_particles;

         // ... else check current size against a (previously) propagated result.
         //     To protect against transient mass exchanges, we only consider situations
         //     where non-dominants are satellites or dominants are centrals
         else if(flag_most_massive==flag_dominant)
            this_subgroup->n_particles_peak=MAX(this_subgroup->n_particles_peak,this_subgroup->n_particles);

         // Propagate peak particle counts forward (unless this halo is a merging fragment)
         tree_horizontal_extended_info *this_subgroup_desc=set_extended_descendant(subgroups,this_subgroup,i_read,n_wrap);
         if(this_subgroup_desc!=NULL){
            if(check_if_type_is_fragmented(this_subgroup->type)==check_if_type_is_fragmented(this_subgroup_desc->type))
               this_subgroup_desc->n_particles_peak=MAX(this_subgroup_desc->n_particles_peak,this_subgroup->n_particles_peak);
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

