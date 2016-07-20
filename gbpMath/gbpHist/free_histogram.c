#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

void free_histogram(hist_info *hist){
  if(check_mode_for_flag(hist->mode,GBP_HISTOGRAM_FIXED)){
     SID_free(SID_FARG hist->bin_count);
  }
  else if(check_mode_for_flag(hist->mode,GBP_HISTOGRAM_IRREGULAR_XLO_DEFINED)){
     SID_free(SID_FARG hist->bin_count);
     SID_free(SID_FARG hist->x_lo);
     SID_free(SID_FARG hist->x_hi);
  }
  else
     SID_trap_error("Invalid mode (%d) specified in free_histogram().",ERROR_LOGIC,hist->mode);
}

