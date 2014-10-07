#include <gbpTrees_build.h>

int set_back_match_file(back_match_info *back_match){
  if(back_match!=NULL){
     if(back_match->halo!=NULL)
        return(back_match->file);
     else
        return(-1);
  }
  else
     return(-1);
}

