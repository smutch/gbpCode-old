#include <gbpTrees_build.h>

tree_horizontal_extended_info *set_extended_descendant(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int i_file,int n_wrap){
  if(halo!=NULL){
     int offset=halo->descendant_file_offset;
     int file  =i_file+offset;
     int index =halo->descendant_index;
     if(offset>0 && index>=0){
        tree_horizontal_extended_info *node_return=&(halos[file%n_wrap][index]);
        if(check_mode_for_flag(node_return->type,TREE_CASE_INVALID))
           SID_trap_error("A descendant pointer points to an invalid halo (%d->%d;%d;%d).",ERROR_LOGIC,i_file,file,index,node_return->type);
        return(node_return);
     }
  }
  return(NULL);
}

