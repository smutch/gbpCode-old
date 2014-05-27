#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_descendant_is_back_matched(tree_horizontal_info *halo,
                                        tree_horizontal_info *descendant){

   //int flag_is_back_matched=FALSE;
   //if(descendant!=NULL){
   //   int              n_back_matches=descendant->n_back_matches;
   //   back_match_info *back_matches  =descendant->back_matches;
   //   for(int i_back_match=0;i_back_match<n_back_matches && !flag_is_back_matched;i_back_match++){
   //      if((back_matches[i_back_match].halo)==halo) flag_is_back_matched=TRUE;
   //   }
   //}
   //return(flag_is_back_matched);
   return((descendant->bridge_backmatch.halo)==halo);
}

