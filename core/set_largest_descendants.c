#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void set_largest_descendants(tree_horizontal_info *halos_i,
                             int                   n_halos_i){
   SID_log("Setting largest descendants...",SID_LOG_OPEN|SID_LOG_TIMER);
   for(int i_halo=0;i_halo<n_halos_i;i_halo++){
      if(check_mode_for_flag(halos_i[i_halo].type,TREE_CASE_MAIN_PROGENITOR)){
         int largest_descendant=0;
         if(halos_i[i_halo].descendant.halo!=NULL){
            largest_descendant=halos_i[i_halo].descendant.halo->n_particles_largest_descendant;
            halos_i[i_halo].n_particles_largest_descendant=MAX(largest_descendant,halos_i[i_halo].n_particles);
         }
         else
            halos_i[i_halo].n_particles_largest_descendant=halos_i[i_halo].n_particles;
      }
      else
         halos_i[i_halo].n_particles_largest_descendant=halos_i[i_halo].n_particles;
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

