#include <gbpTrees_build.h>

void set_extended_dominant_info(tree_horizontal_extended_info **halos,
                                tree_horizontal_extended_info  *this_halo,
                                int flag_parent_has_dominant,
                                int i_file,
                                int index,
                                int n_wrap){

   // Set dominant halo flags at tree leaves, if the halo's parent
   //    does not already have a dominant halo identified
   if(check_mode_for_flag (this_halo->type,TREE_CASE_NO_PROGENITORS) &&
      check_mode_for_flag (this_halo->type,TREE_CASE_MOST_MASSIVE)   &&
      !flag_parent_has_dominant){
      this_halo->type|=TREE_CASE_DOMINANT;
   }

   // Propagate dominant halo flags
   tree_horizontal_extended_info *dom_prog=set_extended_dominant_progenitor(halos,this_halo,n_wrap);
   if(dom_prog!=NULL){
      // Propagate the dominance flag unless a halo changes parent 
      //    and doesn't become the most massive structure
      if(check_mode_for_flag(dom_prog->type,TREE_CASE_DOMINANT) &&
         !(dom_prog->parent_id!=this_halo->parent_id && !check_mode_for_flag(this_halo->type,TREE_CASE_MOST_MASSIVE))){
         this_halo->type|=TREE_CASE_DOMINANT;
      }
   }

   // Set this halo's peak size 
   set_n_particles_peak(this_halo->type,this_halo->n_particles,&(this_halo->n_particles_peak));

   // Set descendant's dominant progenitor.  Descendant inherits this halo's peak size if it is the dominant progenitor
   if(set_extended_new_dominant_progenitor(halos,this_halo,i_file,index,n_wrap)){
      tree_horizontal_extended_info *this_halo_desc=set_extended_descendant(halos,this_halo,i_file,n_wrap);
      this_halo_desc->n_particles_peak=this_halo->n_particles_peak;
   }
}

