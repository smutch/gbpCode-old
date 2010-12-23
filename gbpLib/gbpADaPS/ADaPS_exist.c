#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <gbpCommon.h>
#include <gbpADaPS.h>

int ADaPS_exist(ADaPS *list,
                char  *name_in, ...){
  ADaPS   *current;
  va_list  vargs;
  char     name[ADaPS_NAME_LENGTH];
  va_start(vargs,name_in);
  vsprintf(name,name_in,vargs);
  current=list;
  while(current!=NULL){
    if(!strcmp(name,current->name))
      return(TRUE);
    current=current->next;
  }
  va_end(vargs);
  return(FALSE);
}

