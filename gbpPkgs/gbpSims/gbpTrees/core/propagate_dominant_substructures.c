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
      // Set some group information
      tree_horizontal_extended_info *this_group=&(groups[i_read%n_wrap][i_group]);
      int n_subgroups_group_i=n_subgroups_group[i_read%n_wrap][i_group];
      // Set some defaults
      tree_horizontal_extended_info *dominant_desc=NULL;
      int flag_none_propagated=TRUE;
      int dominant_file_offset=0;
      // Don't bother searching this group if it has no substructures
      //   or if its dominant is dropped for this snapshot
      if(n_subgroups_group_i>0){
         // Initialize: if this group needs to have it's dominant substructure initialized, 
         //             then turn-on the dominant flag for the most massive substructure (MMS)
         tree_horizontal_extended_info *subgroup_MMS=&(subgroups[i_read%n_wrap][i_subgroup]);
         if(check_mode_for_flag(this_group->type,TREE_CASE_NO_PROGENITORS) || 
            check_mode_for_flag(this_group->type,TREE_CASE_REINIT_DOMINANT))
            subgroup_MMS->type|=TREE_CASE_DOMINANT;
         // Propagate: Loop over all of the group's substructures, looking for a dominant descendant
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
                     int this_group_id=this_group->id;
                     if((subgroup_current_desc->parent_id)==this_group_id){
                        dominant_desc       =subgroup_current_desc;
                        dominant_file_offset=subgroup_current->descendant_file_offset;
                     }
                  }
               }
            }
         } // loop over j_subgroup
      }
      // Propagate dominant flag here if we have identified a descendant to propagate it to
      tree_horizontal_extended_info *this_group_desc=set_extended_descendant(groups,this_group,i_read,n_wrap);
      if(dominant_desc!=NULL){
         dominant_desc->type |=TREE_CASE_DOMINANT;
         flag_none_propagated =FALSE;
         // If the dominant skips any snapshots, we have to flag the group as such
         //    so we know not to call for reinitialization in the next snap.
         if(dominant_file_offset>1 && this_group_desc!=NULL)
            this_group_desc->type|=TREE_CASE_DOMINANT_DROPPED;
      }
      // If a dominant flag was not propagated and this group has a descendant, turn-on the re-initialize flag 
      //    or propagate the dropped dominant flag.  Only set a reinitialize flag if the group's descendant is in this
      //    group's main progenitor line.
      if(this_group_desc!=NULL && flag_none_propagated){
         if(check_mode_for_flag(this_group->type,TREE_CASE_DOMINANT_DROPPED))
            this_group_desc->type|=TREE_CASE_DOMINANT_DROPPED;
         else if(this_group->id==this_group_desc->id)
            this_group_desc->type|=TREE_CASE_REINIT_DOMINANT;
      }
      // If a halo was propagated, then the dropped dominant flag should be turned off
      if(!flag_none_propagated)
         this_group->type&=(~TREE_CASE_DOMINANT_DROPPED);
      // Turn-off all re-initialize flags.  All groups have been processed now and we don't want it to be present anywhere in the final output.
      this_group->type&=(~TREE_CASE_REINIT_DOMINANT);
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

