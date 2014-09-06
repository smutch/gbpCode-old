#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

void finalize_histogram(hist_info *hist){
  SID_Allreduce(SID_IN_PLACE, (hist->bin_count), hist->n_bins,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE,&(hist->count_hist),1,           SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE,&(hist->count_all), 1,           SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
}

