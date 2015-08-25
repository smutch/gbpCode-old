#include <gbpTrees_build.h>

tree_horizontal_extended_info *set_extended_dominant_progenitor(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int n_wrap){
  if(halo!=NULL){
     int file  =halo->dominant_progenitor_file;
     int index =halo->dominant_progenitor_index;
     if(file>=0 && index>=0)
        return(&(halos[file%n_wrap][index]));
  }
  return(NULL);
}

