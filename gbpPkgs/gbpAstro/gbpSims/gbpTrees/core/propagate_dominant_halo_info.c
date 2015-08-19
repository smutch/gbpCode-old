#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_dominant_halo_info(tree_horizontal_extended_info **groups,   int *n_groups,
                                  tree_horizontal_extended_info **subgroups,int *n_subgroups,
                                  int        **n_subgroups_group,
                                  int          i_read, // tree snapshot index
                                  int          j_read, // actual snapshot index
                                  int          l_read,
                                  int          i_read_step,
                                  int          n_wrap){
   SID_log("Propagating dominant halo info for snapshot #%03d...",SID_LOG_OPEN,j_read);
   // Process groups
   int i_group;
   int i_subgroup;
   int j_subgroup;
   int flag_returned;
   for(i_group=0,i_subgroup=0;i_group<n_groups[l_read];i_group++){
      tree_horizontal_extended_info *this_group     =&(groups[i_read%n_wrap][i_group]);
      tree_horizontal_extended_info *this_group_desc=set_extended_descendant(groups,this_group,i_read,n_wrap);

      // Set dominant halo flag and peak halo size
      if(check_mode_for_flag(this_group->type,TREE_CASE_NO_PROGENITORS) && !check_mode_for_flag(this_group->type,TREE_CASE_FRAGMENTED_NEW)){
         this_group->type            |=TREE_CASE_DOMINANT;
         this_group->n_particles_peak =this_group->n_particles;
      }

      // Propagate ...
      if(this_group_desc!=NULL){
         // ... dominant halo flags 
         if(check_mode_for_flag(this_group->type,TREE_CASE_DOMINANT))
            this_group_desc->type|=TREE_CASE_DOMINANT;
         // ... n_particles_peak
         this_group_desc->n_particles_peak=MAX(this_group->n_particles_peak,this_group_desc->n_particles);
      }

      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         tree_horizontal_extended_info *this_subgroup     =&(subgroups[i_read%n_wrap][i_subgroup]);
         tree_horizontal_extended_info *this_subgroup_desc=set_extended_descendant(subgroups,this_subgroup,i_read,n_wrap);
         // Set dominant halo flag and peak halo size
         if(check_mode_for_flag(this_subgroup->type,TREE_CASE_NO_PROGENITORS) && !check_mode_for_flag(this_subgroup->type,TREE_CASE_FRAGMENTED_NEW)){
            this_subgroup->type            |=TREE_CASE_DOMINANT;
            this_subgroup->n_particles_peak =this_subgroup->n_particles;
         }

         // Propagate ...
         if(this_subgroup_desc!=NULL){
            // ... dominant halo flags 
            if(check_mode_for_flag(this_subgroup->type,TREE_CASE_DOMINANT)){
               // Don't propagate TREE_CASE_DOMINANT flag if the halo changes parent 
               //    and doesn't become the most massive structure
               if(!(this_subgroup->parent_id!=this_subgroup_desc->parent_id && (this_subgroup_desc->substructure_index!=0)))
                  this_subgroup_desc->type|=TREE_CASE_DOMINANT;
            }
            // ... n_particles_peak
            if(check_mode_for_flag(this_subgroup->type,TREE_CASE_DOMINANT))
               this_subgroup_desc->n_particles_peak=MAX(this_subgroup->n_particles_peak,this_subgroup_desc->n_particles);
            else if((this_subgroup->substructure_index!=0) && (this_subgroup_desc->substructure_index==0))
               this_subgroup_desc->n_particles_peak=this_subgroup->n_particles_peak;
            else
               this_subgroup_desc->n_particles_peak=MAX(this_subgroup->n_particles_peak,this_subgroup_desc->n_particles);
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

