#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

void add_to_histogram(hist_info *hist,double x){
  if(hist->flag_finalized)
     SID_trap_error("An addition has been attempted on an already-finalized histogram in add_to_histogram().",ERROR_LOGIC);
  int index=calc_histogram_index(hist,x);
  if(is_histogram_index_in_range(hist,index)){
     hist->bin_count[index]++;
     hist->count_hist++;
  }
  hist->count_all++;
}
 
