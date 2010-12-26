#include <stdio.h>
#include <common.h>
void free_dtable(dtable *table){
  dremark *next;
  while(table->columns!=NULL){
    remove_dcolumn(table,table->columns->name);
  }
  while(table->remarks!=NULL){
    next=table->remarks->next;
    free(table->remarks->text);
    free(table->remarks);
    table->remarks=next;
  }
  table->n_data=0;
}
