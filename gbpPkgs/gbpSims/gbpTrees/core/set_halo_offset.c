#include <gbpTrees_build.h>

int set_halo_offset(tree_horizontal_info *halo,tree_horizontal_info *target_halo){
  if(halo!=NULL && target_halo!=NULL){
     return((target_halo->file)-(halo->file));
  }
  else
     return(-1);
}

