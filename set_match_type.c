#include <gbpTrees.h>

int set_match_type(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
        return(match->halo->type);
     else
        return(-1);
  }
  else
     return(-1);
}

