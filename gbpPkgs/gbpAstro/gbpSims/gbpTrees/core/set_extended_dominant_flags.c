#include <gbpTrees_build.h>

void set_extended_dominant_flags(tree_horizontal_extended_info *this_halo,tree_horizontal_extended_info *this_halo_desc){
   // Set dominant halo flags
   if(check_mode_for_flag(this_halo->type,TREE_CASE_NO_PROGENITORS) && !check_mode_for_flag(this_halo->type,TREE_CASE_FRAGMENTED_NEW))
      this_halo->type|=TREE_CASE_DOMINANT;
   if(this_halo_desc!=NULL){
      if(check_mode_for_flag(this_halo->type,TREE_CASE_DOMINANT)){
         // Don't propagate TREE_CASE_DOMINANT flag if the halo changes parent 
         //    and doesn't become the most massive structure
         if(!(this_halo->parent_id!=this_halo_desc->parent_id && (this_halo_desc->substructure_index!=0)))
            this_halo_desc->type|=TREE_CASE_DOMINANT;
      }
   }
}
