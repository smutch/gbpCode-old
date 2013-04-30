#include <gbpTrees.h>

int set_match_index(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
        return(match->halo->index);
     else
        return(-1);
  }
  else
     return(-1);
}

