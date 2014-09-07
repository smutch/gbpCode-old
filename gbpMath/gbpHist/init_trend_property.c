#include <stdio.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHist.h>

void init_trend_property(trend_property_info **property,
                         const char           *name,
                         int                   n_hist,
                         void                 *params,
                         void (*init_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         void (*free_function)(trend_property_info *property,void *params_init,int i_hist,int *mode,gbp_va_list *vargs_gbp),
                         int  (*calc_function)(trend_property_info *property,hist_info *hist,void *params_calc)){
   (*property)=(trend_property_info *)SID_malloc(sizeof(trend_property_info));
   strcpy((*property)->name,name);
   (*property)->params       =params;
   (*property)->init_function=init_function;
   (*property)->free_function=free_function;
   (*property)->calc_function=calc_function;
   (*property)->next         =NULL;

   // Create the histogram(s)
   (*property)->n_hist=n_hist;
   (*property)->hist  =(hist_info *)SID_malloc(sizeof(hist_info)*n_hist);
   for(int i_hist=0;i_hist<n_hist;i_hist++){
      int         mode_hist;
      gbp_va_list vargs_gbp;
      init_function((*property),params,i_hist,&mode_hist,&vargs_gbp);
      init_histogram(&((*property)->hist[i_hist]),mode_hist|GBP_HISTOGRAM_GBP_VARGS,&vargs_gbp);
      free_function((*property),params,i_hist,&mode_hist,&vargs_gbp);
   }
}

