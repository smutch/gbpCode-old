#include <gbpTrees_build.h>

int set_match_file(match_info *match){
  if(match!=NULL){
     if(match->halo!=NULL)
        return(match->halo->file);
     else
        return(-1);
  }
  else
     return(-1);
}

