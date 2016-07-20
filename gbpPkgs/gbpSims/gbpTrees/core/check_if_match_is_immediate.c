#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

// Note that the halo pointer in the match_info pointer sent to this
//   function must point to the halos matched FROM, not the target
//   halo matched TO.
int check_if_match_is_immediate(tree_horizontal_info *target_halo,
                                match_info           *match){
   int flag_valid=FALSE;
   if(match!=NULL){
      if(match->halo!=NULL){
         int file_offset=(target_halo->file)-(match->halo->file);
         if(file_offset==1 || file_offset==-1)
            flag_valid=TRUE;
      }
   }
   return(flag_valid);
}
