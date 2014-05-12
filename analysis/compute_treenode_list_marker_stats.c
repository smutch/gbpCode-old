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
#include <gbpTrees_analysis.h>

void compute_treenode_list_marker_stats(tree_info *trees,treenode_list_info *list,tree_markers_stats_info *stats){

   int *hist=(int *)SID_calloc(sizeof(int)*trees->n_snaps);
   for(int i_list=0;i_list<list->n_list;i_list++){
      tree_markers_info markers;
      find_treenode_markers(trees,list->list[i_list],&markers);
      if(markers.half_peak_mass!=NULL){
         hist[markers.half_peak_mass->snap_tree]++;
      }
   }

   int     n_bins=trees->n_snaps;
   int     sum=0;
   for(int i_bin=0;i_bin<n_bins;i_bin++)
      sum+=hist[i_bin];
   size_t *hist_index=NULL;
   merge_sort(hist,(size_t)n_bins,&hist_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
   int i_peak =hist_index[n_bins-1];
   int i_68_lo=hist_index[n_bins-1];
   int i_68_hi=hist_index[n_bins-1];
   int target =(int)(0.68*(double)sum);
   int accum  =0;
   int i_bin  =0;
   while(accum<=target && i_bin<n_bins){
      size_t idx_i=hist_index[n_bins-i_bin-1];
      if(idx_i<i_68_lo) i_68_lo=idx_i;
      if(idx_i>i_68_hi) i_68_hi=idx_i;
      accum+=hist[idx_i];
      i_bin++;
   }
   int i_95_lo=i_68_lo;
   int i_95_hi=i_68_hi;
   target=(int)(0.95*(double)sum);
   while(accum<=target && i_bin<n_bins){
      size_t idx_i=hist_index[n_bins-i_bin-1];
      if(idx_i<i_95_lo) i_95_lo=idx_i;
      if(idx_i>i_95_hi) i_95_hi=idx_i;
      accum+=hist[idx_i];
      i_bin++;
   }
   SID_free(SID_FARG hist);
   SID_free(SID_FARG hist_index);

   stats->t_half_peak_mass_ranges[0][0]=trees->t_list[i_68_lo];
   stats->t_half_peak_mass_ranges[0][1]=trees->t_list[i_68_hi];
   stats->t_half_peak_mass_ranges[1][0]=trees->t_list[i_95_lo];
   stats->t_half_peak_mass_ranges[1][1]=trees->t_list[i_95_hi];
   stats->t_half_peak_mass_peak        =trees->t_list[i_peak];

}

