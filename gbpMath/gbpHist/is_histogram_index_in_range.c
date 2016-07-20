#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

int is_histogram_index_in_range(hist_info *hist,int index){
  return((index>=0) && (index<hist->n_bins));
}

