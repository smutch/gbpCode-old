#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpHist.h>

void finalize_trend(trend_info *trend){
   finalize_trend_property(trend->ordinate);
   trend_property_info *current_coordinate=trend->coordinate_first;
   while(current_coordinate!=NULL){
      finalize_trend_property(current_coordinate);
      current_coordinate=current_coordinate->next;
   }
}

