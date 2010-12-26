#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>
#include <ADaM.h>

void read_dtable_2part(char   *filename_defn,
		       char   *filename_data,
		       dtable *table){

  FILE *fp_defn;
  FILE *fp_data;
  int   n_columns;
  int   n_data;
  int   i,j;
  char *line;
  char  column_name[256];
  char  column_type_txt[256];
  int   column_type;
  int   line_length=0;
  char *data;

  init_dtable(table);

  /* Open files */
  fp_defn=fopen(filename_defn,"r");
  fp_data=fopen(filename_data,"r");

  /* Count lines */
  n_columns=count_lines(fp_defn);
  n_data   =count_lines_data(fp_data);

  /* Grab each column line */
  for(i=0;i<n_columns;i++){
    grab_next_line(fp_defn,&line,&line_length);

    /* Set column name */
    grab_char(line,1,column_name);

    /* Set column type */
    if(count_words(line)>1){
      grab_char(line,2,column_type_txt);
      if(!strcmp(column_type_txt,"int")){
	column_type=ADaM_INT;
	data       =(char *)malloc(sizeof(int)*n_data);
      }
      else if(!strcmp(column_type_txt,"long")){
	column_type=ADaM_LONG;
	data       =(char *)malloc(sizeof(long)*n_data);
      }
      else if(!strcmp(column_type_txt,"double")){
	column_type=ADaM_DOUBLE;
	data       =(char *)malloc(sizeof(double)*n_data);
      }
      else if(!strcmp(column_type_txt,"float")){
	column_type=ADaM_FLOAT;
	data       =(char *)malloc(sizeof(float)*n_data);
      }
      else
	fprintf(stderr,"WARNING: unknown column type {%s}\n",column_type_txt);
    }
    else{
      column_type=ADaM_DOUBLE;
      data       =(char *)malloc(sizeof(double)*n_data);
    }

    /* Read data */
    rewind(fp_data);
    for(j=0;j<n_data;j++){
      grab_next_line_data(fp_data,&line,&line_length);
      switch(column_type){
      case ADaM_INT:
	grab_int(line,i+1,&(((int *)data)[j]));
	break;
      case ADaM_LONG:
	grab_long(line,i+1,&(((long *)data)[j]));
	break;
      case ADaM_FLOAT:
	grab_float(line,i+1,&(((float *)data)[j]));
	break;
      case ADaM_DOUBLE:
	grab_double(line,i+1,&(((double *)data)[j])); 
	break;
      }
    }

    /* Store column */
    add_dcolumn(table,
		column_name,
		column_type,
		(void *)data,
		n_data);

    free(data);
  }

  fclose(fp_defn);
  fclose(fp_data);
  free(line);
}
