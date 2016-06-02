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
            // By making a first check on file_root, we effectively
            //    avoid turning the halo with the longest descendant
            //    line into a fragmented halo.
            if(old_match->halo->file_root>new_match->halo->file_root)
               flag_valid=FALSE;
            else if(old_match->halo->file_root==new_match->halo->file_root){
               if(old_match->halo->file_root>=new_match->halo->file_root)
                  flag_valid=FALSE;
            }
         }
      }
   }
   return(flag_valid);
}
