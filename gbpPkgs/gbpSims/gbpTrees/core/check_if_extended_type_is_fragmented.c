#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_extended_type_is_fragmented(tree_horizontal_extended_info *halo){
   return(check_mode_for_flag(halo->type,TREE_CASE_FRAGMENTED_STRAYED)||
          check_mode_for_flag(halo->type,TREE_CASE_FRAGMENTED_NORMAL) ||
          check_mode_for_flag(halo->type,TREE_CASE_FRAGMENTED_EJECTED));
}

