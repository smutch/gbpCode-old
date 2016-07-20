#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_halo_is_better_main_progenitor(tree_horizontal_info *target_halo,
                                            match_info           *old_progenitor,
                                            match_info           *new_progenitor){
   // We tried using n_particles to set the main progenitor but this does not lead to the proper behavior
   //    during complicated many-body mergers.  Core-satellite switching events lead to mistakes in the
   //    assignment of centrals to the tree currently occupying the central core, leading to a long-lived 
   //    structure assembled from segments of dropped bits that should belong to other trees.
   //int file_offset_new=(target_halo->file)-(new_progenitor->halo->file);
   //int file_offset_old=(target_halo->file)-(old_progenitor->halo->file);
   //return(file_offset_new<file_offset_old || (file_offset_new==file_offset_old && new_progenitor->halo->n_particles>old_progenitor->halo->n_particles));

   // Use matching score to decide the best main progenitor.
   //int file_offset_new=(target_halo->file)-(new_progenitor->halo->file);
   //int file_offset_old=(target_halo->file)-(old_progenitor->halo->file);
   //return(file_offset_new<file_offset_old || (file_offset_new==file_offset_old && (new_progenitor->score)>(old_progenitor->score)));

   // Use best match criterion to decide the best main progenitor.  This does 
   //   a slightly better job in some cases due to the 2-way match check.
   if(check_if_match_is_better(target_halo,old_progenitor,new_progenitor));
}

