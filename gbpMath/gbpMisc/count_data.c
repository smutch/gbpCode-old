#include <stdio.h>
#include <string.h>
#include <common.h>
int count_data(FILE *fp){
  int   n_data=0;
  char *line  =NULL;
  int   n     =0;
  int   i;
  int   r;
  int   flag_data;
  int   flag_check;
  while(!feof(fp)){
    r=getline(&line,&n,fp);
    if(r>0) {
      if
      flag_data =FALSE;
      flag_check=TRUE;
      i=0; while(i<strlen(line) && (temp_char)
      if(flag_data)
	n_data++;
      free(line);
      n=0;
    }
  }
  rewind(fp);
  return(n_data);
}
