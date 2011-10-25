#include <stdio.h>
#include <gbpADaPS.h>

void ADaPS_free(void **list){
  ADaPS *current;
  ADaPS *next;
  current=(*list);
  while(current!=NULL){
    next=current->next;
    ADaPS_deallocate(&current);
    current=next;
  }
  (*list)=NULL;
}
