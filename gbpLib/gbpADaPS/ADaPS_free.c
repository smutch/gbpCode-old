#include <stdio.h>
#include <gbpADaPS.h>

void ADaPS_free(void **list){
  ADaPS *current;
  ADaPS *next;
  current=(ADaPS *)(*list);
  while(current!=NULL){
    next=current->next;
    #if USE_DEBUGGER 
       SID_log("Freeing {%s}...",SID_LOG_OPEN,current->name);
    #endif
    ADaPS_deallocate(&current);
    #if USE_DEBUGGER
       SID_log("Done.",SID_LOG_CLOSE);
    #endif
    current=next;
  }
  (*list)=NULL;
}
