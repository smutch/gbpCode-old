#include <gbpLib.h>
#include <gbpTrees_build.h>

int check_validity_of_tree_case_flag(int flag){
   return(flag>=0 && flag!=TREE_CASE_INVALID && flag!=TREE_CASE_UNPROCESSED);
}

