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
      tree_horizontal_extended_info *dominant_desc=NULL;
      // Set some group information
      tree_horizontal_extended_info *this_group=&(groups[i_read%n_wrap][i_group]);
      int n_subgroups_group_i =n_subgroups_group[i_read%n_wrap][i_group];
      int flag_none_propagated=TRUE;
      // Don't bother if there are no substructures in the group
      if(n_subgroups_group_i>0){
         // Initialize: if this group needs to have it's dominant substructure initialized, 
         //             then turn-on the dominant flag for the most massive substructure (MMS)
         tree_horizontal_extended_info *subgroup_MMS=&(subgroups[i_read%n_wrap][i_subgroup]);
         if(check_mode_for_flag(this_group->type,TREE_CASE_NO_PROGENITORS) || 
            check_mode_for_flag(this_group->type,TREE_CASE_REINIT_DOMINANT)){
            subgroup_MMS->type  |=TREE_CASE_DOMINANT;
            flag_none_propagated =TRUE; // We have initialized a group.  If it isn't propagated, we'll have to do so again in its next snapshot.
         }
         // Propagate: Loop over all of the group's substructures, looking for a dominant descendant
         int this_group_id=this_group->id;
         for(int j_subgroup=0;j_subgroup<n_subgroups_group_i;i_subgroup++,j_subgroup++){
            // Remembering that the subgroups are ordered by current size, stop when we
            //    find the first dominant descendant.  Keep looping to keep counters straight.
            if(dominant_desc==NULL){
               tree_horizontal_extended_info *subgroup_current=&(subgroups[i_read%n_wrap][i_subgroup]);
               // If this subgroup is marked dominant ...
               if(check_mode_for_flag(subgroup_current->type,TREE_CASE_DOMINANT)){
                  // ... and it has a valid descendant ...
                  tree_horizontal_extended_info *subgroup_current_desc=set_extended_descendant(subgroups,subgroup_current,i_read,n_wrap);
                  if(subgroup_current_desc!=NULL){
                     // ... and that descendant is in the group's main progenitor line ...
                     if((subgroup_current_desc->parent_id)==this_group_id)
                        dominant_desc=subgroup_current_desc;
                  }
               }
            }
         } // loop over j_subgroup
      }
      // Any groups that were supposed to be initialized here will have to be done in the group's next snapshot
      else if(check_mode_for_flag(this_group->type,TREE_CASE_NO_PROGENITORS) || 
              check_mode_for_flag(this_group->type,TREE_CASE_REINIT_DOMINANT))
         flag_none_propagated=TRUE;
      // Propagate dominant flag here.
      if(dominant_desc!=NULL){
         dominant_desc->type |=TREE_CASE_DOMINANT;
         flag_none_propagated =FALSE;
      }
      // Turn-on the re-initialize flag for the group descendant if we failed to propagate a dominant substructure or if
      //    this halo was marked as needing to be initialized but doesn't have any substructures.  Do so though,
      //    only if the descendant group is part of this group's main progenitor line.
      tree_horizontal_extended_info *this_group_desc=set_extended_descendant(groups,this_group,i_read,n_wrap);
      if(this_group_desc!=NULL && flag_none_propagated){
         if(this_group->id==this_group_desc->id)
            this_group_desc->type|=TREE_CASE_REINIT_DOMINANT;
      }
      // Turn-off all re-initialize flags.  All groups have been processed now and we don't want it to be present anywhere in the final output.
      this_group->type&=(~TREE_CASE_REINIT_DOMINANT);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

