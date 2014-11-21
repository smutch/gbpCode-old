#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_validity_of_main_progenitor(match_info *old_MP,
                                      match_info *new_MP){
   int flag_valid=TRUE;
   // If the old halo is a 2-way match ...
   if(check_mode_for_flag(old_MP->halo->type,TREE_CASE_2WAY_MATCH)){
      // ... and the new one is not, don't use it.
      if(!check_mode_for_flag(new_MP->halo->type,TREE_CASE_2WAY_MATCH))
         flag_valid=FALSE;
      // ... else, keep the one with the largest number of particles
      else if((old_MP->halo->n_particles)>(new_MP->halo->n_particles))
         flag_valid=FALSE;
      else
         flag_valid=TRUE;
   }
   // ... else if the new halo is a 2-way match, use it.
   else if(check_mode_for_flag(new_MP->halo->type,TREE_CASE_2WAY_MATCH))
      flag_valid=TRUE;
   // ... else, keep the one with the largest number of particles
   else if((old_MP->halo->n_particles)>(new_MP->halo->n_particles))
      flag_valid=FALSE;
   else
      flag_valid=TRUE;
   return(flag_valid);
}

