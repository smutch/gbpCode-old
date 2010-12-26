#include <stdio.h>
#include <string.h>
#include <common.h>
int exist_dcolumn(dtable *table,
		  char    name[DCOLUMN_NAME_LENGTH]){
  dcolumn *current;
  current=table->columns;
  while(current!=NULL){
    if(!strcmp(name,current->name))
      return(TRUE);
    current=current->next;
  }
  return(FALSE);
}

