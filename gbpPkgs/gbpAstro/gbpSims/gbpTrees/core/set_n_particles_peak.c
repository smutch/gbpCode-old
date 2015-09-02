#include <gbpTrees_build.h>

void set_n_particles_peak(int type,int n_particles_halo,int *n_particles_peak){
   // Set peak halo size 
   if(check_mode_for_flag(type,TREE_CASE_NO_PROGENITORS))
      (*n_particles_peak)=n_particles_halo;
   // ... increment inherited value unless a most massive substructure is not dominant
   else if(!(!check_mode_for_flag(type,TREE_CASE_DOMINANT) && check_mode_for_flag(type,TREE_CASE_MOST_MASSIVE)))
      (*n_particles_peak)=MAX((*n_particles_peak),n_particles_halo);
}

