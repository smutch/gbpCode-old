#include <gbpTrees_build.h>

void set_extended_n_particles_peak(tree_horizontal_extended_info *this_halo,tree_horizontal_extended_info *this_halo_desc){
   // Set peak halo size
   if(check_mode_for_flag(this_halo->type,TREE_CASE_NO_PROGENITORS))
      this_halo->n_particles_peak =this_halo->n_particles;

   // Propagate peak size ...
   if(this_halo_desc!=NULL){
      // ... if the halo is dominant, then accumulate peak
      if(check_mode_for_flag(this_halo->type,TREE_CASE_DOMINANT))
         this_halo_desc->n_particles_peak=MAX(this_halo->n_particles_peak,MAX(this_halo_desc->n_particles_peak,this_halo_desc->n_particles));
      // ... if not, and the descendant becomes the most massive substructure (MMS), then freeze peak
      else if(this_halo_desc->substructure_index==0)
         this_halo_desc->n_particles_peak=MAX(this_halo->n_particles_peak,this_halo_desc->n_particles_peak);
      // ... if not dominant and descendant is not MMS, then accumulate peak
      else
         this_halo_desc->n_particles_peak=MAX(this_halo->n_particles_peak,MAX(this_halo_desc->n_particles_peak,this_halo_desc->n_particles));
   }
}
