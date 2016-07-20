#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

double histogram_bin_x_mid(hist_info *hist,int bin){
  return(0.5*(histogram_bin_x_lo(hist,bin)+histogram_bin_x_hi(hist,bin)));
}

