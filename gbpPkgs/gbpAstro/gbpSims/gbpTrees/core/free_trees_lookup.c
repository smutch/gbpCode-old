#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void free_trees_lookup(tree_info *trees){
   int i_wrap;
   if(trees->group_indices!=NULL){
      for(i_wrap=0;i_wrap<trees->n_wrap_lookup;i_wrap++)
         SID_free(SID_FARG trees->group_indices[i_wrap]);
      SID_free(SID_FARG trees->group_indices);
   }
   if(trees->group_array!=NULL){
      for(i_wrap=0;i_wrap<trees->n_wrap_lookup;i_wrap++)
         SID_free(SID_FARG trees->group_array[i_wrap]);
      SID_free(SID_FARG trees->group_array);
   }
   if(trees->subgroup_indices!=NULL){
      for(i_wrap=0;i_wrap<trees->n_wrap_lookup;i_wrap++)
         SID_free(SID_FARG trees->subgroup_indices[i_wrap]);
      SID_free(SID_FARG trees->subgroup_indices);
   }
   if(trees->subgroup_array!=NULL){
      for(i_wrap=0;i_wrap<trees->n_wrap_lookup;i_wrap++)
         SID_free(SID_FARG trees->subgroup_array[i_wrap]);
      SID_free(SID_FARG trees->subgroup_array);
   }
}

