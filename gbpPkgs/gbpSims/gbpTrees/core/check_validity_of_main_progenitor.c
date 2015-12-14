#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_validity_of_main_progenitor(tree_horizontal_info *descendant_halo,
                                      match_info           *old_MP,
                                      match_info           *new_MP){
   int flag_valid=TRUE;
   if(old_MP!=NULL){
      if(old_MP->halo!=NULL){
         if(new_MP==NULL)
            flag_valid=FALSE;
         else if(new_MP->halo==NULL)
            flag_valid=FALSE;
         else{
            int file_offset_new=(descendant_halo->file)-(new_MP->halo->file);
            int file_offset_old=(descendant_halo->file)-(old_MP->halo->file);
            // Always keep main progenitors as close to the descendant as possible
            if(file_offset_new>file_offset_old)
               flag_valid=FALSE;
            // We use n_particles here instead of score, because otherwise small merging halos
            //    can be tagged as the main progenitor if they sink efficiently to the centre
            //    of a large system and the matching is weighted strongly to the centre.
            else if((old_MP->halo->n_particles)>(new_MP->halo->n_particles))
               flag_valid=FALSE;
         }
      }
   }
   return(flag_valid);
}
