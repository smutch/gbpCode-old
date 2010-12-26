#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ADaM.h>
#include <common.h>
int add_dremark(dtable  *table,
		char    *text){
  dremark *last;
  dremark *next;
  dremark *new_dremark;

  new_dremark      =(dremark *)malloc(sizeof(dremark));
  new_dremark->text=(char    *)malloc(sizeof(char)*(strlen(text)+1));
  strcpy(new_dremark->text,text);
  new_dremark->next=NULL;
  next=table->remarks;
  last=table->remarks;
  while(next!=NULL){
    last=next;
    next=last->next;
  }
  if(last!=next)
    last->next=new_dremark;
  else
    table->remarks=new_dremark;

  return(ERROR_NONE);
}

