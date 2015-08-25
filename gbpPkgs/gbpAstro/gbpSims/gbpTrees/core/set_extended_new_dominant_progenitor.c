#include <gbpTrees_build.h>

int set_extended_new_dominant_progenitor(tree_horizontal_extended_info **halos,
                                         tree_horizontal_extended_info  *halo,
                                         int i_file,
                                         int index,
                                         int n_wrap){
  int flag_new_prog=FALSE;
  if(halo!=NULL){
     tree_horizontal_extended_info *desc=set_extended_descendant(halos,halo,i_file,n_wrap);
     if(desc!=NULL){
        tree_horizontal_extended_info *old_dom=set_extended_dominant_progenitor(halos,desc,n_wrap);
        if(old_dom==NULL)
           flag_new_prog=TRUE;
        // Check here if we have a better dominant progenitor
        else{
           int flag_dom_old=check_mode_for_flag(old_dom->type,TREE_CASE_DOMINANT);
           int flag_dom_new=check_mode_for_flag(halo->type,   TREE_CASE_DOMINANT);
           int n_p_old     =old_dom->n_particles_peak;
           int n_p_new     =halo->n_particles_peak;
           if(flag_dom_old){
              if(flag_dom_new && n_p_new>n_p_old) flag_new_prog=TRUE;
           }
           else{
              if(flag_dom_new || n_p_new>n_p_old) flag_new_prog=TRUE;
           }
        }
     }
     if(flag_new_prog){
        desc->dominant_progenitor_file =i_file;
        desc->dominant_progenitor_index=index;
     }
  }
  return(flag_new_prog);
}

