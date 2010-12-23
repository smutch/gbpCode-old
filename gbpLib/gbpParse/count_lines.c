/********************************/
/*         count_lines.c        */
/*         -------------        */
/* Count the number of lines in */
/* a file pointed to by *fp     */
/*------------------------------*/
/* Author:    Greg Poole        */
/* Last edit: May 29/06         */
/********************************/
#include <stdio.h>
#include <gbpCommon.h>
int count_lines(FILE *fp){
  int   n_lines=0;
  char *line   =NULL;
  int   n      =0;
  int   r;
  while(!feof(fp)){
    r=getline(&line,&n,fp);
    if(r>0) n_lines++;
  }
  rewind(fp);
  return(n_lines);
}
