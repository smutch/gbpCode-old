#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void free_trend(trend_info **trend){
   // Free the ordinate property
   free_trend_property(&((*trend)->ordinate));

   // Free the coordinate properties
   trend_property_info *current_coordinate=(*trend)->coordinate_first;
   while(current_coordinate!=NULL){
      trend_property_info *current_coordinate_next=current_coordinate->next;
      free_trend_property(&current_coordinate);
      current_coordinate=current_coordinate_next;
   }

   SID_free(SID_FARG (*trend));
}

