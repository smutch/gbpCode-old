#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMisc.h>
#include <gbpSort.h>

void compute_histogram_range(hist_info *hist,double confidence_percent,int mode,double *x_peak,double *x_lo,double *x_hi){

   // Sanity check
   if(check_mode_for_flag(mode,GBP_HISTOGRAM_RANGE_ALL) && check_mode_for_flag(mode,GBP_HISTOGRAM_RANGE_HIST))
      SID_trap_error("Multiple range mode flags have been set in compute_histogram_range().",ERROR_LOGIC);
   if(!check_mode_for_flag(mode,GBP_HISTOGRAM_RANGE_ALL) && !check_mode_for_flag(mode,GBP_HISTOGRAM_RANGE_HIST))
      SID_trap_error("No range mode flags have been set in compute_histogram_range().",ERROR_LOGIC);

   // Set defaults
   int i_peak=hist->n_bins/2;
   int i_lo  =0;
   int i_hi  =hist->n_bins-1;

   // Set the population to nomalize by
   size_t n_norm;
   if(check_mode_for_flag(mode,GBP_HISTOGRAM_RANGE_ALL))
      n_norm=hist->count_all;
   else if(check_mode_for_flag(mode,GBP_HISTOGRAM_RANGE_HIST))
      n_norm=hist->count_hist;

   if(n_norm>0){
      // Sort counts
      size_t *hist_index=NULL;
      merge_sort(hist->bin_count,(size_t)(hist->n_bins),&hist_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);

      // Set confidence interval count
      size_t target=(size_t)(0.01*confidence_percent*(double)n_norm);

      // Find range
      size_t accum=0;
      int    i_bin=0;
      i_peak =hist_index[hist->n_bins-1];
      i_lo   =hist_index[hist->n_bins-1];
      i_hi   =hist_index[hist->n_bins-1];
      while(accum<target && i_bin<hist->n_bins){
         size_t idx_i=hist_index[hist->n_bins-i_bin-1];
         if(idx_i<i_lo) i_lo=idx_i;
         if(idx_i>i_hi) i_hi=idx_i;
         accum+=hist->bin_count[idx_i];
         i_bin++;
      }

      // Clean-up
      SID_free(SID_FARG hist_index);
   }

   // Invert indices if need-be
   if(hist->flag_bin_order_inverted){
      int i_temp=i_lo;
      i_lo=i_hi;
      i_hi=i_temp;
   }

   // Convert array indices into values
   (*x_peak)=histogram_bin_x_mid(hist,i_peak);
   (*x_lo)  =histogram_bin_x_lo (hist,i_lo);
   (*x_hi)  =histogram_bin_x_hi (hist,i_hi);

}
 
