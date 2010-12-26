#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>

void init_table_column(dtable *table,
		       char   *column_name,
		       int     column_type){
  dcolumn *new;
  dcolumn *next;
  dcolumn *last;
  dcolumn *current;

  /* Create new column */
  new=(dcolumn *)malloc(sizeof(dcolumn));
  strcpy(new->name,column,name);
  new->type  =column_type;
  new->n_data=0;
  new->data  =NULL;
  new->next  =NULL;

  /* Remove it if it already exists */
  dcolumn_remove(table,column_name);

  /* Add it to the linked list */
  current=table->columns;
  next=table->columns;
  last=table->columns;
  while(next!=NULL){
    last=next;
    next=last->next;
  }
  if(last!=next)
    last->next=new;
  else
    (*list)=new;
}
