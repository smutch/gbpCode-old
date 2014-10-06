#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

void finalize_trend_property(trend_property_info *property){
  for(int i_hist=0;i_hist<property->n_hist;i_hist++) finalize_histogram(&(property->hist[i_hist]));
}

