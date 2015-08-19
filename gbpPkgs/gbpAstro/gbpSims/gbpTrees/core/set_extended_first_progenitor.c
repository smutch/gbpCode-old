#include <gbpTrees_build.h>

tree_horizontal_extended_info *set_extended_first_progenitor(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int n_wrap){
  if(halo!=NULL){
     int file =halo->first_progenitor_file;
     int index=halo->first_progenitor_index;
     if(file>=0 && index>=0)
        return(&(halos[file%n_wrap][index]));
  }
  return(NULL);
}

