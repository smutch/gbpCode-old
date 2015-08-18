#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_halo_is_merger(int type){
   return(!check_mode_for_flag(type,TREE_CASE_MAIN_PROGENITOR) && !check_mode_for_flag(type,TREE_CASE_INVALID));
}

