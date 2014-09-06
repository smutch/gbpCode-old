#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <gbpLib.h>
#include <gbpMisc.h>

double histogram_bin_x_lo(hist_info *hist,int bin){
  double r_val;

  // Sanity check
  if(!is_histogram_index_in_range(hist,bin))
     SID_trap_error("Invalid bin requested (%d) in histogram_bin_x_hi().",ERROR_LOGIC,bin);

  if(check_mode_for_flag(hist->mode,GBP_HISTOGRAM_FIXED)){
     switch(hist->flag_bin_order_inverted){
        case FALSE:
           r_val=hist->x_min+((double)(bin))*hist->dx;
           break;
        case TRUE:
           r_val=hist->x_min+((double)(hist->n_bins-bin-1))*hist->dx;
           break;
     }
  }
  else if(check_mode_for_flag(hist->mode,GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED))
     r_val=hist->x_lo[bin];
  else
     SID_trap_error("Invalid mode (%d) specified in histogram_bin_x_lo().",ERROR_LOGIC,hist->mode);
  return(r_val);
}

