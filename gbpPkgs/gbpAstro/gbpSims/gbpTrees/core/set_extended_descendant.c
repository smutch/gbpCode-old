#include <gbpTrees_build.h>

tree_horizontal_extended_info *set_extended_descendant(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int i_file,int n_wrap){
  if(halo!=NULL){
     int offset=halo->descendant_file_offset;
     int file  =i_file+offset;
     int index =halo->descendant_index;
     if(offset>0 && index>=0)
        return(&(halos[file%n_wrap][index]));
  }
  return(NULL);
}

