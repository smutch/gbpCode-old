#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHist.h>

void init_trend_ordinate(trend_info  *trend,
                         const char  *name,
                         void        *params,
                         void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc)){
   init_trend_property(&(trend->ordinate),name,TRUE,1,params,init_function,free_function,calc_function);
}

