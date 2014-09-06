#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMisc.h>

int calc_histogram_index(hist_info *hist,double x){
  int r_val=-1;
  if(check_mode_for_flag(hist->mode,GBP_HISTOGRAM_FIXED)){
     r_val=(int)((x-hist->x_min)/hist->dx);
     // Deal with bin-order swapping
     switch(hist->flag_bin_order_inverted){
        case TRUE:
           r_val=(hist->n_bins-r_val-1);
           break;
     }
  }
  else if(check_mode_for_flag(hist->mode,GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED)){
     int flag_done=FALSE;
     switch(hist->flag_bin_order_inverted){
        case FALSE:
           for(int i_bin=0;i_bin<hist->n_bins && !flag_done;i_bin++){
              if(x>=histogram_bin_x_lo(hist,i_bin)){
                 if(x<histogram_bin_x_hi(hist,i_bin)){
                    r_val    =i_bin;
                    flag_done=TRUE;
                 }
              }
              else
                 flag_done=TRUE;
           }
           break;
        case TRUE:
           for(int i_bin=(hist->n_bins-1);i_bin>=0 && !flag_done;i_bin++){
              if(x>=histogram_bin_x_lo(hist,i_bin)){
                 if(x<histogram_bin_x_hi(hist,i_bin)){
                    r_val    =i_bin;
                    flag_done=TRUE;
                 }
              }
              else
                 flag_done=TRUE;
           }
           break;
     }
  }
  else
     SID_trap_error("Invalid mode (%d) specified in clear_histogram().",ERROR_LOGIC,hist->mode);

  return(r_val);
}

