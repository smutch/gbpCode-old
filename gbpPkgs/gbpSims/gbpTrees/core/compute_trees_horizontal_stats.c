#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void compute_trees_horizontal_stats(void *halos_in,int n_halos,int n_halos_max,tree_horizontal_stats_info *stats,int flag_write_cases){
   int i_halo;

   // Initialize counters etc
   init_trees_horizontal_stats(stats,n_halos);

   switch(flag_write_cases){
      case TRUE:{
         tree_horizontal_info *halos;
         halos=(tree_horizontal_info *)halos_in;
         for(i_halo=0;i_halo<n_halos_max;i_halo++)
            add_to_trees_horizontal_stats(stats,&(halos[i_halo]));
         }
         break;
      case FALSE:{
         tree_horizontal_extended_info *halos;
         halos=(tree_horizontal_extended_info *)halos_in;
         for(i_halo=0;i_halo<n_halos_max;i_halo++)
            add_to_trees_horizontal_stats(stats,&(halos[i_halo]));
         }
         break;
   }

   if(n_halos!=(n_halos_max-stats->n_invalid))
      SID_trap_error("There is an incorrect number of out-of-bounds halos (i.e. %d!=%d)",ERROR_LOGIC,n_halos,(n_halos_max-stats->n_invalid));
   if(stats->n_unprocessed!=0)
      SID_trap_error("A number of halos (%d) are still marked as unprocessed",ERROR_LOGIC,stats->n_unprocessed);
}

