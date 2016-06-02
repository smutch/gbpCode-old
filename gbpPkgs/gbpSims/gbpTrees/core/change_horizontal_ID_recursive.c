#include <gbpTrees_build.h>

void change_horizontal_ID_recursive(tree_horizontal_info *halo,int i_file,int id_1,int id_2){
   tree_horizontal_info *current;

   // Change IDs and root file numbers here
   if(halo->id==id_1){
      halo->id       =id_2;
      halo->file_root=i_file;
   }

   // Walk the tree
   current=halo->first_progenitor.halo;
   while(current!=NULL){
      change_horizontal_ID_recursive(current,i_file,id_1,id_2);
      current=current->next_progenitor.halo;
   }
}

