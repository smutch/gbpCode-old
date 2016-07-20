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
   for(i_group=0,i_subgroup=0;i_group<n_groups[l_read];i_group++){
      // Set merger flags...
      //    ... first, find the primary halo ...
      tree_horizontal_extended_info *this_group   =&(groups[i_read%n_wrap][i_group]);
      tree_horizontal_extended_info *current_group=set_extended_first_progenitor(groups,this_group,n_wrap);
      tree_horizontal_extended_info *primary_group=current_group;
      int n_progenitors=0;
      while(current_group!=NULL){
         // Decide if this halo is a better primary
         if((current_group->n_particles_peak)>(primary_group->n_particles_peak)){
            // Don't set fragmented halos to be primaries if it can be avoided
            int flag_frag_old=check_if_extended_type_is_fragmented(primary_group);
            int flag_frag_new=check_if_extended_type_is_fragmented(current_group);
            if(!flag_frag_new || flag_frag_new==flag_frag_old)
               primary_group=current_group;
         }
         current_group=set_extended_next_progenitor(groups,current_group,n_wrap);
         n_progenitors++;
      }
      //    ... then, set the appopriate merger flags.
      if(n_progenitors>1){
         current_group=set_extended_first_progenitor(groups,this_group,n_wrap);
         while(current_group!=NULL){
            if(current_group==primary_group)
               current_group->type|=TREE_CASE_MERGER_PRIMARY;
            else 
               current_group->type|=TREE_CASE_MERGER;
            current_group=set_extended_next_progenitor(groups,current_group,n_wrap);
         }
      }
      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group];i_subgroup++,j_subgroup++){
         // Set merger flags...
         //    ... first, find the primary halo ...
         tree_horizontal_extended_info *this_subgroup   =&(subgroups[i_read%n_wrap][i_subgroup]);
         tree_horizontal_extended_info *current_subgroup=set_extended_first_progenitor(subgroups,this_subgroup,n_wrap);
         tree_horizontal_extended_info *primary_subgroup=current_subgroup;
         n_progenitors=0;
         while(current_subgroup!=NULL){
            // Decide if this halo is a better primary
            if((current_subgroup->n_particles_peak)>(primary_subgroup->n_particles_peak)){
               // Don't set fragmented halos to be primaries if it can be avoided
               int flag_frag_old=check_if_extended_type_is_fragmented(primary_subgroup);
               int flag_frag_new=check_if_extended_type_is_fragmented(current_subgroup);
               if(!flag_frag_new || flag_frag_new==flag_frag_old)
                  primary_subgroup=current_subgroup;
            }
            current_subgroup=set_extended_next_progenitor(subgroups,current_subgroup,n_wrap);
            n_progenitors++;
         }
         //    ... then, set the appropriate merger flags.
         if(n_progenitors>1){
            current_subgroup=set_extended_first_progenitor(subgroups,this_subgroup,n_wrap);
            while(current_subgroup!=NULL){
               if(current_subgroup==primary_subgroup)
                  current_subgroup->type|=TREE_CASE_MERGER_PRIMARY;
               else
                  current_subgroup->type|=TREE_CASE_MERGER;
               current_subgroup=set_extended_next_progenitor(subgroups,current_subgroup,n_wrap);
            }
         }
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

