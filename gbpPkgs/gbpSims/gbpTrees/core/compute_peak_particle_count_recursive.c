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

void compute_peak_particle_count_recursive(tree_info *trees,tree_node_info *this_halo,int *flag_fragmented_return,int *n_particles_peak,int *n_particles_inclusive_peak){

  // Initialize some flags
  int flag_most_massive=check_mode_for_flag(this_halo->tree_case,TREE_CASE_MOST_MASSIVE);
  int flag_dominant    =check_mode_for_flag(this_halo->tree_case,TREE_CASE_DOMINANT);
  int flag_fragmented  =check_if_type_is_fragmented(this_halo->tree_case);

  // Initialize result to zero
  this_halo->n_particles_peak          =0;
  this_halo->n_particles_inclusive_peak=0;

  // Loop over each progenitor (recursively)
  tree_node_info *first_progenitor  =this_halo->progenitor_first;
  tree_node_info *current_progenitor=first_progenitor;
  while(current_progenitor!=NULL){
     int flag_fragmented_prog;
     int n_particles_peak_prog;
     int n_particles_inclusive_peak_prog;
     compute_peak_particle_count_recursive(trees,current_progenitor,&flag_fragmented_prog,&n_particles_peak_prog,&n_particles_inclusive_peak_prog);
     // Don't propagate peak particle counts if this is a merging fragment
     if(flag_fragmented==flag_fragmented_prog){
        this_halo->n_particles_peak          =MAX(this_halo->n_particles_peak,          n_particles_peak_prog);
        this_halo->n_particles_inclusive_peak=MAX(this_halo->n_particles_inclusive_peak,n_particles_inclusive_peak_prog);
     }
     current_progenitor=current_progenitor->progenitor_next;
  }

  // Initialize values at leaves
  if(first_progenitor==NULL || this_halo->n_particles_peak==0){
     this_halo->n_particles_peak          =this_halo->n_particles;
     this_halo->n_particles_inclusive_peak=this_halo->n_particles_inclusive;
  }
  // ... else, apply dominant halo logic
  else if(this_halo->parent_top==NULL || flag_most_massive==flag_dominant){
     this_halo->n_particles_peak          =MAX(this_halo->n_particles_peak,          this_halo->n_particles);
     this_halo->n_particles_inclusive_peak=MAX(this_halo->n_particles_inclusive_peak,this_halo->n_particles_inclusive);
  }

  // Pass results up the tree
  if(n_particles_peak!=NULL)
     (*n_particles_peak)=this_halo->n_particles_peak;
  if(n_particles_inclusive_peak!=NULL)
     (*n_particles_inclusive_peak)=this_halo->n_particles_inclusive_peak;
  if(flag_fragmented_return!=NULL)
     (*flag_fragmented_return)=flag_fragmented;
}

