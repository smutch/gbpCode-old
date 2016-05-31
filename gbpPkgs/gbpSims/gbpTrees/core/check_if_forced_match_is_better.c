#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_forced_match_is_better(tree_horizontal_info *target_halo,
                                    match_info           *old_match,
                                    match_info           *new_match){
   int flag_valid=TRUE;
   if(old_match!=NULL){
      if(old_match->halo!=NULL){
         if(new_match==NULL)
            flag_valid=FALSE;
         else if(new_match->halo==NULL)
            flag_valid=FALSE;
         else{
            // First, prioritize 2way matches
            if(old_match->flag_two_way)
               flag_valid=FALSE;
            // ... then, prioritize larger progenitors
            else if(!new_match->flag_two_way && old_match->halo->n_particles>new_match->halo->n_particles)
               flag_valid=FALSE;
         }
      }
   }
   return(flag_valid);
}
