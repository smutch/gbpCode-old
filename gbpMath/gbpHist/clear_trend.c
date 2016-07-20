#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void clear_trend(trend_info **trend){
   // Free the ordinate property
   clear_histogram((*trend)->ordinate->hist);

   // Free the coordinate properties
   trend_property_info *current_coordinate=(*trend)->coordinate_first;
   while(current_coordinate!=NULL){
      trend_property_info *current_coordinate_next=current_coordinate->next;
      for(int i_hist=0;i_hist<current_coordinate->n_hist;i_hist++)
         clear_histogram(&(current_coordinate->hist[i_hist]));
      current_coordinate=current_coordinate_next;
   }
}

