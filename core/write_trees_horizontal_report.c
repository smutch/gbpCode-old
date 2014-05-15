#include <gbpLib.h>
#include <gbpTrees_build.h>

void write_trees_horizontal_report(int                   n_halos_i,
                                   int                   n_halos_max,
                                   tree_horizontal_info *halos_i){
   tree_horizontal_stats_info stats;
   compute_trees_horizontal_stats(halos_i,n_halos_i,n_halos_max,&stats,TRUE);
   SID_log("Results (estimates which may change with continued processing):",SID_LOG_OPEN);
   SID_log("# of halos               = %d",SID_LOG_COMMENT,stats.n_halos);
   SID_log("# of mergers             = %d",SID_LOG_COMMENT,stats.n_mergers);
   SID_log("# of strayed halos       = %d",SID_LOG_COMMENT,stats.n_strayed);
   SID_log("# of dropped halos       = %d",SID_LOG_COMMENT,stats.n_dropped);
   SID_log("# of bridged halos       = %d",SID_LOG_COMMENT,stats.n_bridged);
   SID_log("# of emerged progenitors = %d",SID_LOG_COMMENT,stats.n_emerged_progenitors);
   SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
}

