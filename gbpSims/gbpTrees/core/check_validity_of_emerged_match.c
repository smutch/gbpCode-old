#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_validity_of_emerged_match(tree_horizontal_info *halo_i,
                                    back_match_info      *back_match,
                                    int n_search){
   // Perform some checks to see if we want to make a match to this emerged halo
   int flag_valid=TRUE;

   // 1) Because we may be recursively finding emerged halos, make sure we haven't
   //    exceeded the search interval.  If so, invalidate the match.
   int total_offset=back_match->file-halo_i->file;
   if(total_offset<=0)
      SID_trap_error("Invalid emerged halo snapshot offset (ie. %d<=0).",ERROR_LOGIC,total_offset);
   if(total_offset>n_search)
      flag_valid=FALSE;

   // 2) Check that we are not matching to a descendant of the initial bridged match ...
   //check_if_halo_is_descendant(halo_i,back_match->halo,n_search);

   return(flag_valid);
}

