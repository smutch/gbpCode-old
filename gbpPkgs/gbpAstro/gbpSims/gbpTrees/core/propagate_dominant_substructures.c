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
      // Don't bother if there are no substructures in the group
      if(n_subgroups_group[i_read%n_wrap][i_group]>0){
         // Find the group's dominant progenitor (if there is one)
         int i_subgroup_dominant=i_subgroup; // default to the first/most-massive 
         int flag_found         =FALSE;
         for(int j_subgroup=0;j_subgroup<n_subgroups_group[i_read%n_wrap][i_group] && !flag_found;i_subgroup++,j_subgroup++){
            // Loop over all of the substructure's progenitors looking for halos with dominant flags set
            tree_horizontal_extended_info *current_subgroup     =&(subgroups[i_read%n_wrap][i_subgroup]);
            tree_horizontal_extended_info *current_subgroup_prog=set_extended_first_progenitor(subgroups,current_subgroup,n_wrap);
            while(current_subgroup_prog!=NULL && !flag_found){
               if(check_mode_for_flag(current_subgroup_prog->type,TREE_CASE_DOMINANT)){
                  i_subgroup_dominant=i_subgroup;
                  flag_found         =TRUE; // This way, if two dominant halos come in, we choose the most massive
               }
               current_subgroup_prog=set_extended_next_progenitor(subgroups,current_subgroup_prog,n_wrap);
            }
         }
         // We've chosen a dominant substructure.  Add the flag.
         tree_horizontal_extended_info *subgroup_dominant=&(subgroups[i_read%n_wrap][i_subgroup_dominant]);
         subgroup_dominant->type|=TREE_CASE_DOMINANT;
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

