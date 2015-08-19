#include <gbpTrees_build.h>

tree_horizontal_extended_info *set_extended_descendant(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int i_file,int n_wrap){
  if(halo!=NULL){
     int file =i_file+halo->file_offset;
     int index=halo->index;
     if(file>=0 && index>=0)
        return(&(halos[file%n_wrap][index]));
  }
  return(NULL);
}

