#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_halo_is_fragmented(int type){
   return(check_mode_for_flag(type,TREE_CASE_FRAGMENTED_STRAYED)  ||
          check_mode_for_flag(type,TREE_CASE_FRAGMENTED_RETURNED) ||
          check_mode_for_flag(type,TREE_CASE_FRAGMENTED_EXCHANGED));
}

