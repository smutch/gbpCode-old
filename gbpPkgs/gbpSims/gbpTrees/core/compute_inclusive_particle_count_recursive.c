#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void compute_inclusive_particle_count_recursive(tree_info *trees,tree_node_info *this_halo,int *n_particles_substructure){
  // Loop over each substructure (recursively)
  int n_particles_substructure_sum    =0;
  tree_node_info *first_substructure  =this_halo->substructure_first;
  tree_node_info *current_substructure=first_substructure;
  while(current_substructure!=NULL){
     int n_particles_substructure_i;
     compute_inclusive_particle_count_recursive(trees,current_substructure,&n_particles_substructure_i);
     n_particles_substructure_sum+=n_particles_substructure_i;
     current_substructure         =current_substructure->substructure_next;
  }
  // If this pointer is NULL, then it means that
  //    this function has been passed a group and
  //    we don't want to alter it's particle count,
  //    which is already inclusive by construction
  if(n_particles_substructure!=NULL){
     this_halo->n_particles_inclusive+=n_particles_substructure_sum;
     (*n_particles_substructure)=this_halo->n_particles_inclusive;
  }
}

