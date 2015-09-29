#include <gbpTrees_build.h>

tree_horizontal_extended_info *set_extended_first_progenitor(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int n_wrap){
  if(halo!=NULL){
     int file =halo->first_progenitor_file;
     int index=halo->first_progenitor_index;
     if(file>=0 && index>=0){
        tree_horizontal_extended_info *node_return=&(halos[file%n_wrap][index]);
        if(check_mode_for_flag(node_return->type,TREE_CASE_INVALID))
           SID_trap_error("A first_progenitor pointer points to an invalid halo (->%d;%d;%d).",ERROR_LOGIC,file,index,node_return->type);
        return(node_return);
     }
  }
  return(NULL);
}

