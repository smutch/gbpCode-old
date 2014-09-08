#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void init_trend(trend_info **trend,
                const char  *name,
                void        *params,
                void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc)){
   (*trend)                  =(trend_info *)SID_malloc(sizeof(trend_info));
   (*trend)->ordinate        =NULL;
   (*trend)->coordinate_first=NULL;
   (*trend)->coordinate_last =NULL;

   // Create the ordinate property
   init_trend_ordinate((*trend),name,params,init_function,free_function,calc_function);
}

