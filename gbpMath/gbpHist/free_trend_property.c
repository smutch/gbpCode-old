#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void free_trend_property(trend_property_info **property){
   if((*property)->hist!=NULL){
      for(int i_hist=0;i_hist<(*property)->n_hist;i_hist++)
         free_histogram(&((*property)->hist[i_hist])); 
      SID_free(SID_FARG (*property)->hist);
   }
   SID_free(SID_FARG (*property));
}     

