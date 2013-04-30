#include <gbpTrees.h>

int set_match_snapshot(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
        return(match->halo->snap);
     else
        return(-1);
  }
  else
     return(-1);
}

