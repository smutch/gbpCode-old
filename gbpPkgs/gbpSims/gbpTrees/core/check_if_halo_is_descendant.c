#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

int check_if_halo_is_descendant(tree_horizontal_info *possible_progenitor,
                                tree_horizontal_info *possible_descendant,
                                int n_search){
   int flag_is_progenitor=FALSE;
   int k_file;
   int l_file;
   tree_horizontal_info *current;
   current=possible_progenitor->descendant.halo;
   if(current!=NULL)
      k_file =current->file;
   l_file=k_file;

   // Loop over the main progenitor line of the original bridge match
   while(current!=NULL &&
         k_file>=l_file && k_file<=(possible_descendant->file) && // not true when we reach past the rolling array bounds
         !flag_is_progenitor){
      // ... if we come to the halo we have been asked to match to, invalidate the match.
      if(current==possible_descendant)
        flag_is_progenitor=TRUE;
      // ... if not, keep moving ...
      current=current->descendant.halo;
      l_file=k_file;
      if(current!=NULL)
         k_file=current->file;
   }
   return(flag_is_progenitor);
}

