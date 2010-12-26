#include <stdio.h>
#include <string.h>
#include <common.h>
int remove_dcolumn(dtable *table, 
		   char    name[DCOLUMN_NAME_LENGTH]){
  dcolumn *last   =NULL;
  dcolumn *current=NULL;
  dcolumn *remove =NULL;
  current=table->columns;
  while(current!=NULL && remove==NULL){
    if(!strcmp(name,current->name)){
      if(last==NULL){
        remove =table->columns;
        table->columns=table->columns->next;        
      }
      else{
        remove    =current;
        last->next=current->next;
      }
    }
    last   =current;
    current=last->next;
  }
  if(remove!=NULL){
    free_dcolumn(&remove);
    table->n_columns--;
  }
  return(ERROR_NONE);
}
