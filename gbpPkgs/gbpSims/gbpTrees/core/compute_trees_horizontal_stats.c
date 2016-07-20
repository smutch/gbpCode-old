#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void compute_trees_horizontal_stats(void *halos_in,int n_halos,int n_halos_max,tree_horizontal_stats_info *stats,int flag_write_cases){
   // Initialize counters etc
   init_trees_horizontal_stats(stats,n_halos);

   switch(flag_write_cases){
      case TRUE:{
         tree_horizontal_info *halos=(tree_horizontal_info *)halos_in;
         for(int i_halo=0;i_halo<n_halos_max;i_halo++){
            tree_horizontal_info *halo_i=&(halos[i_halo]);
            int id=halo_i->id;
            int descendant_id=-1;
            if(halo_i->descendant.halo!=NULL)
               descendant_id=halo_i->descendant.halo->id;
            int n_particles=halo_i->n_particles;
            int type       =halo_i->type;
            add_to_trees_horizontal_stats(stats,id,descendant_id,n_particles,type);
         }
         }
         break;
      case FALSE:{
         tree_horizontal_extended_info *halos=(tree_horizontal_extended_info *)halos_in;
         for(int i_halo=0;i_halo<n_halos_max;i_halo++){
            tree_horizontal_extended_info *halo_i=&(halos[i_halo]);
            int id           =halo_i->id;
            int descendant_id=halo_i->descendant_id;
            int n_particles  =halo_i->n_particles;
            int type         =halo_i->type;
            add_to_trees_horizontal_stats(stats,id,descendant_id,n_particles,type);
         }
         }
         break;
   }

   if(n_halos!=(n_halos_max-stats->n_invalid))
      SID_trap_error("There is an incorrect number of out-of-bounds halos (i.e. %d!=%d)",ERROR_LOGIC,n_halos,(n_halos_max-stats->n_invalid));
   if(stats->n_unprocessed!=0)
      SID_trap_error("A number of halos (%d) are still marked as unprocessed",ERROR_LOGIC,stats->n_unprocessed);
}

