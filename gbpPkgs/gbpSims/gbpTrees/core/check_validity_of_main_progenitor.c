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
            // This works well for centrally weighted matches.  Using n_particles instead
            //    has been checked but does not work as well.  Unfortunately, this choice
            //    comes at the expense of sometimes selecting a small substructure, that 
            //    sinks quickly to (or passes transiently through) the middle of a larger 
            //    halo, as the main progenitor.
            else if((old_MP->score)>(new_MP->score))
               flag_valid=FALSE;
         }
      }
   }
   return(flag_valid);
}
