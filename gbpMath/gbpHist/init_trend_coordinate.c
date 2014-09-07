#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHist.h>

void init_trend_coordinate(trend_info  *trend,
                           const char  *name,
                           void        *params,
                           void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                           void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                           int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc)){
   if(trend->ordinate==NULL)
      SID_trap_error("A trend coordinate is being set before its ordinate has been defined.",ERROR_LOGIC);

   // Create new coordinate
   trend_property_info *coordinate_new;
   init_trend_property(&coordinate_new,name,trend->ordinate->hist->n_bins,params,init_function,free_function,calc_function);

   // Add it to the linked list
   if(trend->coordinate_first!=NULL)
      trend->coordinate_last->next=coordinate_new;
   else
      trend->coordinate_first=coordinate_new;
   trend->coordinate_last=coordinate_new;
}

