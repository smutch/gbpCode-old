#include <stdio.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpHist.h>

int add_item_to_trend(trend_info *trend,int mode,void *item_ordinate,...){
   // Parse the mode
   va_list  vargs;
   va_start(vargs,item_ordinate);
   void *item_coordinate;
   if(check_mode_for_flag(mode,GBP_ADD_ITEM_TO_TREND_COORDINATE))
      item_coordinate=item_ordinate;
   else if(check_mode_for_flag(mode,GBP_ADD_ITEM_TO_TREND_ORDINATE))
      item_coordinate=(void *)va_arg(vargs,void *);
   else
      SID_trap_error("Invalid mode (%d) passed to add_item_to_trend.",ERROR_LOGIC,mode);

   int i_ordinate=trend->ordinate->calc_function(trend->ordinate,&(trend->ordinate->hist[0]),item_ordinate);
   if(is_histogram_index_in_range(&(trend->ordinate->hist[0]),i_ordinate)){
      add_to_histogram_index(&(trend->ordinate->hist[0]),i_ordinate);
      // Loop over the coordinates
      trend_property_info *current_coordinate=trend->coordinate_first;
      while(current_coordinate!=NULL){
         hist_info *hist_i      =&(current_coordinate->hist[i_ordinate]);
         int        i_coordinate=current_coordinate->calc_function(current_coordinate,hist_i,item_coordinate);
         add_to_histogram_index(hist_i,i_coordinate);
         current_coordinate=current_coordinate->next;
      }
   }
   else
      i_ordinate=-1;

   va_end(vargs);
   return(i_ordinate);
}

