#include <gbpLib.h>
#include <gbpTrees_build.h>

void add_substructure_info(tree_horizontal_info *halos,
                           int                  *n_subgroups_group,
                           int                  *n_particles_groups,
                           int                   n_groups,
                           int                   n_subgroups,
                           int                   flag_match_halos){
   int i_group   =0;
   int i_subgroup=0;
   for(i_group=0;i_group<n_groups;i_group++){
      if(flag_match_halos==MATCH_SUBGROUPS){
         for(int j_subgroup=0;j_subgroup<n_subgroups_group[i_group];j_subgroup++,i_subgroup++){
            halos[i_subgroup].n_particles_parent =n_particles_groups[i_group];
            if(j_subgroup==0)
               halos[i_subgroup].type|=TREE_CASE_MOST_MASSIVE;
         }
      }
   }
   if(flag_match_halos==MATCH_SUBGROUPS && i_subgroup!=n_subgroups)
      SID_trap_error("Subgroup counts don't match (ie %d!=%d;%d)in add_substructure_info().",ERROR_LOGIC,i_subgroup,n_subgroups,n_groups);
}
