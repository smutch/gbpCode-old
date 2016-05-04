#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_match_is_better(tree_horizontal_info *target_halo,
                             match_info           *old_match,
                             match_info           *new_match){
   int flag_valid=TRUE;
   if(old_match!=NULL){
      if(old_match->halo!=NULL){
         if(new_match==NULL)
            flag_valid=FALSE;
         else if(new_match->halo==NULL)
            flag_valid=FALSE;
         else{
            int file_offset_new=(target_halo->file)-(new_match->halo->file);
            int file_offset_old=(target_halo->file)-(old_match->halo->file);
            // Since we might use this code to check back matches,
            //    correct the sence of the offset
            if(file_offset_new<0){
               file_offset_new=-file_offset_new;
               file_offset_old=-file_offset_old;
            }
            // Sanity check
            if(file_offset_new<=0 || file_offset_old<=0)
               SID_trap_error("There is something odd about the file offsets (%d->%d and %d->%d) in check_if_match_is_better().",ERROR_LOGIC,
                              new_match->halo->file,target_halo->file,old_match->halo->file,target_halo->file);
            // The best match is always the most immediate
            if(file_offset_new>file_offset_old)
               flag_valid=FALSE;
            // ... otherwise, select via secondary criteria ...
            else if(file_offset_new==file_offset_old){
               // The following XXX'd-out scheme fails when snapshot densities are high.  In these cases, small
               //   halo's passing transiently through cores spend several snapshots in this state, creating
               //   2-way matches and leading to the problem this choice was meant to correct.  For this reason,
               //   we now simply take the match with the largest score.
               // XXX // With f_good>0, the following scheme was found to create lots of poor matches due to small
               // XXX //   halos getting matched to/by the centrals of large systems, leading to lots of spurious
               // XXX //   and very large fragmented halos, particularly when snapshot cadence becomes dense.  For
               // XXX //   that reason, f_good must be kept to zero so that all possible back matches created by 
               // XXX //   FoF core swaps are considered.
               // XXX // If equally immediate, give first preference to 2-way matches ...
               // XXX if(old_match->flag_two_way)
               // XXX    flag_valid=FALSE;
               // XXX // ... otherwise, give secondary preference to the halo with the best metrics.
               // XXX //     Unfortunately, the best choice of metric(s) is not perfectly obvious:
               // XXX //    ... using match score for this decision introduces the problem
               // XXX //        of sometimes selecting a small substructure that sinks quickly to (or
               // XXX //        passes transiently through) the middle of a larger halo.  This 
               // XXX //        problem is reduced by giving preference to 2-way matches (above).
               // XXX //    ... using n_particles for this decision can be compromised
               // XXX //        by ambiguities introduced by the large halo size fluctuations
               // XXX //        introduced by transient exchanges of outer halo particles or
               // XXX //        (equivilantly) by changes in the identity of the substructure
               // XXX //        core the halo finder reports as the centre of the system
               // XXX else if(!new_match->flag_two_way && (old_match->score>new_match->score))
               // XXX    flag_valid=FALSE;
               if(old_match->score>new_match->score)
                  flag_valid=FALSE;
            }
         }
      }
   }
   return(flag_valid);
}
