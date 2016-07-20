#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_type_is_fragmented(int type){
   return(check_mode_for_flag(type,TREE_CASE_FRAGMENTED_STRAYED)||
          check_mode_for_flag(type,TREE_CASE_FRAGMENTED_NORMAL) ||
          check_mode_for_flag(type,TREE_CASE_FRAGMENTED_OTHER));
}

