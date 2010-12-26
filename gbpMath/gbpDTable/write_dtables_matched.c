#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>
#include <ADaM.h>

void write_dtables_matched(dtable *table1,
			   dtable *table2,
			   char   *filename,
			   char   *match_name){

  FILE    *fp;
  int      i,j,k;
  dremark *current_remark;
  dcolumn *current_column;
  int     *match_data;
  int      flag_match=FALSE;
  int      match_column=0;
  int      last_column;
  int      flag;

  fp=fopen(filename,"w");

  /* Set-up matching (if needed) */
  if(strcmp(match_name,"") && match_name!=NULL){
    current_column=table1->columns;
    flag=TRUE;
    i=1;
    while(current_column!=NULL){
      if(!strcmp(current_column->name,match_name)){
	if(current_column->type!=ADaM_INT && current_column->type!=ADaM_LONG){
	  fprintf(stderr,"WARNING: invalid match column {type %d instead of %d}...ignoring.\n",
		  current_column->type,
		  ADaM_INT);
	}
	else{
	  match_column=i;
	  match_data  =(int *)current_column->data;
	  flag_match  =TRUE;
	}
      }
      else
	last_column=i;
      i++;
      current_column=current_column->next;
    }
    current_column=table2->columns;
    while(current_column!=NULL){
      last_column=i;
      i++;
      current_column=current_column->next;
    }
    if(!flag_match){
      fprintf(stderr,"WARNING: match column {%s} not found ... dumping all.\n",match_name);
    }
  }
  else
    last_column=table1->n_columns+table2->n_columns;

  /****************/
  /* Write header */
  /****************/

  /* Remarks */
  current_remark=table1->remarks;
  while(current_remark!=NULL){
    fprintf(fp,"# %s\n",current_remark->text);
    current_remark=current_remark->next;
  }

  /* Column info from 1st table*/
  current_column=table1->columns;
  i=1;
  while(current_column!=NULL){
    if(!flag_match || (flag_match && i!=match_column))
      fprintf(fp,"# Column %2d %s\n",i++,current_column->name);
    current_column=current_column->next;
  }

  /* Column info from 2nd table*/
  current_column=table2->columns;
  while(current_column!=NULL){
    fprintf(fp,"# Column %2d %s\n",i++,current_column->name);
    current_column=current_column->next;
  }

  /**************/
  /* Write data */
  /**************/

  /* Loop over all lines of data */
  for(j=0;j<table1->n_data;j++){

    /* Check for matching criteria */
    if(flag_match)
      if(!(match_data[j]))
	continue;

    /* Loop over all columns in 1st table */
    current_column=table1->columns;
    i=1;
    while(current_column!=NULL){
      if(!flag_match || (flag_match && i!=match_column)){
	switch(current_column->type){
	case ADaM_INT:
	  fprintf(fp,"%d",((int *)(current_column->data))[j]);
	  break;
	case ADaM_LONG:
	  fprintf(fp,"%ld",((long *)(current_column->data))[j]);
	  break;
	case ADaM_FLOAT:
	  fprintf(fp,"%e",((float *)(current_column->data))[j]);
	  break;
	case ADaM_DOUBLE:
	  fprintf(fp,"%le",((double *)(current_column->data))[j]);
	  break;
	}
	if(i!=last_column)
	  fprintf(fp,"  ");
      }
      if(i==last_column)
	fprintf(fp,"\n");
      i++;
      current_column=current_column->next;
    }
    
    /* Loop over all columns in 2nd table */
    current_column=table2->columns;
    k=match_data[j];
    while(current_column!=NULL){
      switch(current_column->type){
      case ADaM_INT:
	fprintf(fp,"%d",((int *)(current_column->data))[k]);
	break;
      case ADaM_LONG:
	fprintf(fp,"%ld",((long *)(current_column->data))[k]);
	break;
      case ADaM_FLOAT:
	fprintf(fp,"%e",((float *)(current_column->data))[k]);
	break;
      case ADaM_DOUBLE:
	fprintf(fp,"%le",((double *)(current_column->data))[k]);
	break;
      }
      if(i!=last_column)
	fprintf(fp,"  ");
      else
	fprintf(fp,"\n");
      i++;
      current_column=current_column->next;
    }

  }

  fclose(fp);
}
