#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <gbpCommon.h>
#include <gbpADaPS.h>

void *ADaPS_fetch(ADaPS *list,
                  char  *name_in,...){
  ADaPS   *current;
  int      flag;
  void    *r_val;
  va_list  vargs;
  char     name[ADaPS_NAME_LENGTH];
  va_start(vargs,name_in);
  vsprintf(name,name_in,vargs);
  current=list;
  r_val  =NULL;
  flag   =TRUE;
  while(current!=NULL && flag){
    if(!strcmp(name,current->name)){
      r_val=current->data;
      flag =FALSE;
    }
    current=current->next;
  }
  va_end(vargs);
  return(r_val);
}

