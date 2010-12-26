#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>
#include <ADaM.h>

void read_table_data(char   *filename,
		     dtable *table){

  FILE    *fp;
  int      n_data;
  int      i;
  char    *line;
  char     column_name[256];
  char     column_type_txt[256];
  int      column_type;
  int      line_length=0;
  dcolumn *current;

  /* Open file and count the number of data lines */
  fp    =fopen(filename,"r");
  n_data=count_lines_data(fp);

  /* Allocate data for all columns */
  current=table->columns;
  while(current!=NULL){
    current->n=n_data;
    switch(current->type){
    case ADaM_LONG:
      current->data=(void *)malloc(sizeof(long)*n_data);
      break;
    case ADaM_INT:
      current->data=(void *)malloc(sizeof(int)*n_data);
      break;
    case ADaM_DOUBLE:
      current->data=(void *)malloc(sizeof(double)*n_data);
      break;
    case ADaM_FLOAT:
      current->data=(void *)malloc(sizeof(float)*n_data);
      break;
    }
    current=current->next;
  }

  grab_next_line(fp,&line,&line_length);
  for(i=0;i<n_data;i++){
    while(!check_comment(line))
      grab_next_line(fp,&line,&line_length);
    current=table->columns;
    while(current!=NULL){
      switch(current->type){
      case ADaM_LONG:
	grab_long(line,current->column,&(((long *)current->data)[i_data]));
	break;
      case ADaM_INT:
	grab_int(line,current->column,&(((int *)current->data)[i_data]));
	break;
      case ADaM_DOUBLE:
	grab_float(line,current->column,&(((double *)current->data)[i_data]));
	break;
      case ADaM_FLOAT:
	grab_double(line,current->column,&(((float *)current->data)[i_data]));
	break;
      }
      current=current->next;
    }
  }
  fclose(fp);
}
