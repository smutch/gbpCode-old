/**********************************/
/*         grab_nth_line.c        */
/*         ---------------        */
/* Return the nth line in the     */
/* file pointed to by *fp         */
/*--------------------------------*/
/* Author:    Greg Poole          */
/* Last edit: May 29/06           */
/**********************************/
#ifndef MAX_LINE_LENGTH
#define MAX_LINE_LENGTH 7500
#endif
#ifndef ERROR_NONE
#define ERROR_NONE 0
#endif
#ifndef ERROR_LINE_TOO_LONG
#define ERROR_LINE_TOO_LONG 101
#endif
#ifndef ERROR_FILE_TOO_SHORT
#define ERROR_FILE_TOO_SHORT 102
#endif
#include <stdio.h>
#include <string.h>
int grab_nth_line(FILE *fp,int n,char *line){
  int  counter;
  int  error=ERROR_NONE;
  counter=0;
  while(counter<n){
    if(fscanf(fp,"%s\n",line)!=EOF)
      counter++;
    else{
      counter=n+1;
      error=ERROR_FILE_TOO_SHORT;
    }
  }
  if(error==ERROR_NONE && !(strlen(line)<MAX_LINE_LENGTH))
    error=ERROR_LINE_TOO_LONG;
  return(error);
}
