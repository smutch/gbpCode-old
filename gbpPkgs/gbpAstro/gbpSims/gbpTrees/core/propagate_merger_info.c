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
      int group_id                =groups[i_read%n_wrap][i_group].id;
      int group_n_particles       =groups[i_read%n_wrap][i_group].n_particles;
      int group_tree_id           =groups[i_read%n_wrap][i_group].tree_id;
      int group_descendant_id     =groups[i_read%n_wrap][i_group].descendant_id;
      int group_type              =groups[i_read%n_wrap][i_group].type;
      int group_file_offset       =groups[i_read%n_wrap][i_group].file_offset;
      int group_index             =groups[i_read%n_wrap][i_group].index;
      int group_n_particles_parent=groups[i_read%n_wrap][i_group].n_particles_parent;
      int group_n_particles_desc  =groups[i_read%n_wrap][i_group].n_particles_desc;
      int group_n_particles_proj  =groups[i_read%n_wrap][i_group].n_particles_proj;
      int group_score_desc        =groups[i_read%n_wrap][i_group].score_desc;
      int group_score_prog        =groups[i_read%n_wrap][i_group].score_prog;
      int group_snap_bridge       =groups[i_read%n_wrap][i_group].snap_bridge;
      int group_file_bridge       =groups[i_read%n_wrap][i_group].file_bridge;
      int group_index_bridge      =groups[i_read%n_wrap][i_group].index_bridge;
      int group_id_bridge         =groups[i_read%n_wrap][i_group].id_bridge;


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

