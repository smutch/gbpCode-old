#include <gbpTrees.h>

int set_match_id(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
        return(match->halo->id);
     else
        return(-1);
  }
  else
     return(-1);
}

