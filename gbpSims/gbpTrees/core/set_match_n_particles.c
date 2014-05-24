#include <gbpTrees_build.h>

int set_match_n_particles(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
       return(match->halo->n_particles);
     else
       return(-1);
  }
  else
     return(-1);
}

