#include <stdio.h>
#include <gbpCommon.h>
#include <gbpADaPS.h>

void ADaPS_deallocate(ADaPS **remove){
  SID_free((void **)(&((*remove)->data)));
  SID_free((void **)(remove));
}
