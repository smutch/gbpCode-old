#include <stdio.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpADaPS.h>

void ADaPS_status(ADaPS *list){
  ADaPS *current;
  current=list;
  while(current!=NULL){
    fprintf(stderr,"%30s %2d\n",
            current->name,
            current->mode);
    current=current->next;
  }
}
