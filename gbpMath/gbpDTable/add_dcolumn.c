#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ADaM.h>
#include <common.h>
int add_dcolumn(dtable  *table,
		char    *name,
		int      type,
		void    *data,
		int      n_data){
  dcolumn *last;
  dcolumn *next;
  dcolumn *new_dcolumn;
  int      i;
  if(table->n_data==0 || (table->n_data>0 && table->n_data==n_data)){
    remove_dcolumn(table,name);
    new_dcolumn=(dcolumn *)malloc(sizeof(dcolumn));
    new_dcolumn->n_data=n_data;
    strcpy(new_dcolumn->name,name);
    new_dcolumn->type=type;
    if(type==ADaM_DOUBLE){
      new_dcolumn->data=(void *)malloc(sizeof(double)*n_data);
      for(i=0;i<n_data;i++)
	((double *)new_dcolumn->data)[i]=((double *)data)[i];
    }
    else if(type==ADaM_FLOAT){
      new_dcolumn->data=(void *)malloc(sizeof(float)*n_data);
      for(i=0;i<n_data;i++)
	((float *)new_dcolumn->data)[i]=((float *)data)[i];
    }
    else if(type==ADaM_LONG){
      new_dcolumn->data=(void *)malloc(sizeof(long)*n_data);
      for(i=0;i<n_data;i++)
	((long *)new_dcolumn->data)[i]=((long *)data)[i];
    }
    else if(type==ADaM_INT){
      new_dcolumn->data=(void *)malloc(sizeof(int)*n_data);
      for(i=0;i<n_data;i++)
	((int *)new_dcolumn->data)[i]=((int *)data)[i];
    }
    new_dcolumn->next=NULL;
    next=table->columns;
    last=table->columns;
    while(next!=NULL){
      last=next;
      next=last->next;
    }
    if(last!=next)
      last->next=new_dcolumn;
    else{
      table->columns=new_dcolumn;
      table->n_data =n_data;
    }
    table->n_columns++;
    return(ERROR_NONE);
  }
  else{
    fprintf(stderr,"ERROR: discrepant data column size {%d instead of %d}!\n",
            n_data,
            table->n_data);
    return(-1);
  }
}

