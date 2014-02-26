#include <stdio.h>
#include <gbpCommon.h>
#include <gbpADaPS.h>

void ADaPS_deallocate(ADaPS **remove){
  if((*remove)->free_function!=NULL)
     ((*remove)->free_function)(SID_FARG (*remove)->data,(*remove)->free_function_params);
  else
     SID_free(SID_FARG (*remove)->data);
  SID_free(SID_FARG (*remove));
}

