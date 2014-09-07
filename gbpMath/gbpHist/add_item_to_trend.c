#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void add_item_to_trend(trend_info *trend,void *item){
   int i_ordinate=trend->ordinate->calc_function(trend->ordinate,trend->ordinate->hist,item);
   if(is_histogram_index_in_range(trend->ordinate->hist,i_ordinate)){
      add_to_histogram_index(trend->ordinate->hist,i_ordinate);
      // Loop over the coordinates
      trend_property_info *current_coordinate=trend->coordinate_first;
      while(current_coordinate!=NULL){
         hist_info *hist_i      =&(current_coordinate->hist[i_ordinate]);
         int        i_coordinate=current_coordinate->calc_function(current_coordinate,hist_i,item);
         add_to_histogram_index(hist_i,i_coordinate);
         current_coordinate=current_coordinate->next;
      }
   }
}

