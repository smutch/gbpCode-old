#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

// Note that this function is rather inefficient.  The scan over a descendant's progenitors
//   has a lot of redundant calculations, since the exact same calculation gets repeated by all 
//   its progenitors.  It would be better if this calculation was done once for the 
//   descendants and the information passed down to the progenitors.
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
      tree_horizontal_extended_info *this_group     =&(groups[i_read%n_wrap][i_group]);
      tree_horizontal_extended_info *this_group_desc=set_extended_descendant(groups,this_group,i_read,n_wrap);
      if(this_group_desc!=NULL){
         // Set merger flags...
         //    ... first, find the primary halo ...
         tree_horizontal_extended_info *current_group;
         tree_horizontal_extended_info *primary_group;
         current_group=set_extended_first_progenitor(groups,this_group_desc,n_wrap);
         primary_group=current_group;
         int n_progenitors=0;
         while(current_group!=NULL){
            // Decide if this halo is a better primary
            if((current_group->n_particles_peak)>(primary_group->n_particles_peak))
               primary_group=current_group;
            current_group=set_extended_next_progenitor(groups,current_group,n_wrap);
            n_progenitors++;
         }

         //    ... set the appopriate merger flags.
         if(n_progenitors>1){
            current_group=set_extended_first_progenitor(groups,this_group_desc,n_wrap);
            primary_group=current_group;
            while(current_group!=NULL){
               if(current_group==primary_group)
                  current_group->type|=TREE_CASE_MERGER_PRIMARY;
               else 
                  current_group->type|=TREE_CASE_MERGER;
               current_group=set_extended_next_progenitor(groups,current_group,n_wrap);
            }
         }
      }

      // Process subgroups
      for(j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_subgroup];i_subgroup++,j_subgroup++){
         tree_horizontal_extended_info *this_subgroup     =&(subgroups[i_read%n_wrap][i_subgroup]);
         tree_horizontal_extended_info *this_subgroup_desc=set_extended_descendant(subgroups,this_subgroup,i_read,n_wrap);
         if(this_subgroup_desc!=NULL){
            // Set merger flags...
            //    ... first, find the primary halo ...
            tree_horizontal_extended_info *current_subgroup;
            tree_horizontal_extended_info *primary_subgroup;
            current_subgroup=set_extended_first_progenitor(subgroups,this_subgroup_desc,n_wrap);
            primary_subgroup=current_subgroup;
            int n_progenitors=0;
            while(current_subgroup!=NULL){
               // Decide if this halo is a better primary
               if((current_subgroup->n_particles_peak)>(primary_subgroup->n_particles_peak))
                  primary_subgroup=current_subgroup;
               current_subgroup=set_extended_next_progenitor(subgroups,current_subgroup,n_wrap);
               n_progenitors++;
            }

            //    ... set the appopriate merger flags.
            if(n_progenitors>1){
               current_subgroup=set_extended_first_progenitor(subgroups,this_subgroup_desc,n_wrap);
               primary_subgroup=current_subgroup;
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
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

