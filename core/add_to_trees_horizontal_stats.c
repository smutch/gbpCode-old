#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees_build.h>

void add_to_trees_horizontal_stats(tree_horizontal_stats_info *stats,int id,int type,int n_particles){

   // Find maximum id
   stats->max_id=MAX(stats->max_id,id);

   // Compute statistcs for mergers
   if(check_mode_for_flag(type,TREE_CASE_MERGER))
      stats->n_mergers++;

   // Compute statistcs for strays
   if(check_mode_for_flag(type,TREE_CASE_STRAYED)){
      stats->max_strayed_size=MAX(stats->max_strayed_size,n_particles);
      stats->n_strayed++;
   }

   // Compute statistcs for dropped halos
   if(check_mode_for_flag(type,TREE_CASE_DROPPED)){
      stats->max_dropped_size=MAX(stats->max_dropped_size,n_particles);
      stats->n_dropped++;
   }

   // Compute statistcs for bridged and emerged halos
   if(check_mode_for_flag(type,TREE_CASE_BRIDGED)){
      stats->max_bridged_size=MAX(stats->max_bridged_size,n_particles);
      stats->n_bridged++;
   }
   if(check_mode_for_flag(type,TREE_CASE_MATCHED_TO_BRIDGE)){
      stats->max_bridge_progenitor_size=MAX(stats->max_bridge_progenitor_size,n_particles);
      stats->n_bridge_progenitors++;
      if(!check_mode_for_flag(type,TREE_CASE_BRIDGE_DEFAULT)){
         stats->max_emerged_progenitor_size=MAX(stats->max_emerged_progenitor_size,n_particles);
         stats->n_emerged_progenitors++;
      }
   }
   if(check_mode_for_flag(type,TREE_CASE_EMERGED)){
      stats->max_emerged_size=MAX(stats->max_emerged_size,n_particles);
      stats->n_emerged++;
   }
   if(check_mode_for_flag(type,TREE_CASE_FRAGMENTED_STRAYED)){
      stats->max_fragmented_strayed_size=MAX(stats->max_fragmented_strayed_size,n_particles);
      stats->n_fragmented_strayed++;
   }
   if(check_mode_for_flag(type,TREE_CASE_FRAGMENTED_RETURNED)){
      stats->max_fragmented_returned_size=MAX(stats->max_fragmented_returned_size,n_particles);
      stats->n_fragmented_returned++;
   }
   if(check_mode_for_flag(type,TREE_CASE_FRAGMENTED_EXCHANGED)){
      stats->max_fragmented_exchanged_size=MAX(stats->max_fragmented_exchanged_size,n_particles);
      stats->n_fragmented_exchanged++;
   }

   // Sanity checks
   if(check_mode_for_flag(type,TREE_CASE_INVALID)){
      stats->n_invalid++;
      if((type-TREE_CASE_INVALID)!=0)
         SID_trap_error("An out-of-bounds halo has been manipulated (type=%d)",ERROR_LOGIC,type);
   }
   if(check_mode_for_flag(type,TREE_CASE_UNPROCESSED))
      stats->n_unprocessed++;
}

