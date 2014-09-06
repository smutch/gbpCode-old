#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMisc.h>

void clear_histogram(hist_info *hist){
   for(int i_bin=0;i_bin<hist->n_bins;i_bin++)
      hist->bin_count[i_bin]=0;
   hist->count_hist=0;
   hist->count_all =0;
}

