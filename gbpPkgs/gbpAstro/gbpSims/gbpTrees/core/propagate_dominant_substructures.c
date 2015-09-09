#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void propagate_dominant_substructures(tree_horizontal_extended_info **groups,   int *n_groups,
                                      tree_horizontal_extended_info **subgroups,int *n_subgroups,
                                      int        **n_subgroups_group,
                                      int          i_read, // tree snapshot index
                                      int          j_read, // actual snapshot index
                                      int          l_read,
                                      int          i_read_step,
                                      int          n_wrap){
   SID_log("Propagating dominant substructures for snapshot #%03d...",SID_LOG_OPEN,j_read);
   // Process groups
   int i_subgroup=0;
   for(int i_group=0;i_group<n_groups[l_read];i_group++){
      tree_horizontal_extended_info *this_group     =&(groups[i_read%n_wrap][i_group]);
      tree_horizontal_extended_info *this_group_desc=set_extended_descendant(groups,this_group,i_read,n_wrap);
      int this_group_id=this_group->id;
      // Don't bother if there are no substructures in the group
      int n_subgroups_group_i=n_subgroups_group[i_read%n_wrap][i_group];
      if(n_subgroups_group_i>0){
         tree_horizontal_extended_info *subgroup_dominant=NULL;
         // Loop over all of the group's substructures, looking for a better dominant substructure than the (default) most massive one
         int flag_stop=FALSE; // Turn on this flag if we've found the dominant halo we want and want to stop the search immediately.
         for(int j_subgroup=0;j_subgroup<n_subgroups_group_i && !flag_stop;i_subgroup++,j_subgroup++){
            tree_horizontal_extended_info *subgroup_current=&(subgroups[i_read%n_wrap][i_subgroup]);
            // Ignore fragmented halos
            if(!check_mode_for_flag(subgroup_current->type,TREE_CASE_FRAGMENTED_STRAYED)  &&
               !check_mode_for_flag(subgroup_current->type,TREE_CASE_FRAGMENTED_RETURNED) &&
               !check_mode_for_flag(subgroup_current->type,TREE_CASE_FRAGMENTED_EXCHANGED)){
               // Loop over all of the substructure's progenitors looking for halos with dominant flags set
               tree_horizontal_extended_info *subgroup_current_progenitor=set_extended_first_progenitor(subgroups,subgroup_current,n_wrap);
               int flag_has_dom_prog=FALSE;
               while(subgroup_current_progenitor!=NULL && !flag_has_dom_prog){
                  if(check_mode_for_flag(subgroup_current_progenitor->type,TREE_CASE_DOMINANT))
                     flag_has_dom_prog=TRUE;
                  else
                     subgroup_current_progenitor=set_extended_next_progenitor(subgroups,subgroup_current_progenitor,n_wrap);
               }
               // If this subgroup has a dominant progenitor ...
               if(flag_has_dom_prog){
                  // ... which was a member of this group's main progenitor line, then use it instead of the default.
                  if((subgroup_current_progenitor->parent_id)==this_group_id)
                     subgroup_dominant=&(subgroups[i_read%n_wrap][i_subgroup]);
               }
            }
         }
         // If we've chosen a dominant substructure.  Add the flag.
         if(subgroup_dominant!=NULL)
            subgroup_dominant->type|=TREE_CASE_DOMINANT;
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

