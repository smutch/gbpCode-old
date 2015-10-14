#include <gbpTrees_build.h>
float set_back_match_score(back_match_info *back_match){
  if(back_match!=NULL){
     if(back_match->halo!=NULL)
        return(back_match->score);
     else
        return(0.);
  }
  else
     return(0.);
}

