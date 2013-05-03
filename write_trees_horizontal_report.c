#include <gbpLib.h>
#include <gbpTrees.h>

void write_trees_horizontal_report(int                   n_halos_i,
                                   int                   n_halos_max,
                                   tree_horizontal_info *halos_i){
   tree_horizontal_stats_info stats;
   compute_trees_horizontal_stats(halos_i,n_halos_i,n_halos_max,&stats,TRUE);
   SID_log("Results (estimates which may change with continued processing):",SID_LOG_OPEN);
   SID_log("# of halos                  =%-8d",SID_LOG_COMMENT,stats.n_halos);
   SID_log("# of simple matches         =%-8d (%d mergers)",SID_LOG_COMMENT,stats.n_simple,stats.n_mergers);
   if(stats.n_strayed>0)
      SID_log("# of strayed halos          =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_strayed,stats.max_strayed_size);
   else
      SID_log("# of strayed halos          =%-8d",SID_LOG_COMMENT,stats.n_strayed);
   if(stats.n_sputtered>0)
      SID_log("# of sputtering halos       =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_sputtered,stats.max_sputtered_size);
   else
      SID_log("# of sputtering halos       =%-8d",SID_LOG_COMMENT,stats.n_sputtered);
   if(stats.n_dropped>0)
      SID_log("# of dropped halos          =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_dropped,stats.max_dropped_size);
   else
      SID_log("# of dropped halos          =%-8d",SID_LOG_COMMENT,stats.n_dropped);
   SID_log("# of bridged halos          =%-8d",SID_LOG_COMMENT,stats.n_bridged);
   SID_log("# of bridge progenitors     =%-8d",SID_LOG_COMMENT,stats.n_bridge_progenitors);
   if(stats.n_emerged_progenitors>0)
      SID_log("# of emerged progenitors    =%-8d (largest=%d particles)",SID_LOG_COMMENT,stats.n_emerged_progenitors,stats.max_emerged_progenitor_size);
   else
      SID_log("# of emerged progenitors    =%-8d",SID_LOG_COMMENT,stats.n_emerged_progenitors);
   SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
}

