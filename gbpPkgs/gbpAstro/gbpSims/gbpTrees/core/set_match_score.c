#include <gbpTrees_build.h>
float set_match_score(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
        return(match->score);
     else
        return(0.);
  }
  else
     return(0.);
}

