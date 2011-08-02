#include <stdio.h>
#include <gbpCommon.h>
#include <gbpADaPS.h>

void ADaPS_deallocate(ADaPS **remove){
  ((*remove)->free_function)((void **)(&((*remove)->data)));
  SID_free(SID_FARG (*remove));
}
