#include <gbpTrees_build.h>

void init_trees_horizontal_stats(tree_horizontal_stats_info *stats,int n_halos){
   stats->n_halos                    =n_halos;
   stats->n_mergers                  =0;
   stats->n_strayed                  =0;
   stats->n_dropped                  =0;
   stats->n_bridged                  =0;
   stats->n_emerged                  =0;
   stats->n_fragmented_strayed       =0;
   stats->n_fragmented_normal        =0;
   stats->n_fragmented_ejected       =0;
   stats->n_emerged_progenitors      =0;
   stats->n_invalid                  =0;
   stats->n_unprocessed              =0;
   stats->max_strayed_size           =0;
   stats->max_dropped_size           =0;
   stats->max_bridged_size           =0;
   stats->max_emerged_size           =0;
   stats->max_fragmented_strayed_size=0;
   stats->max_fragmented_normal_size =0;
   stats->max_fragmented_ejected_size=0;
   stats->max_emerged_progenitor_size=0;
   stats->max_id                     =0;
}

