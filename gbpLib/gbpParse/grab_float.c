#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse.h>

int grab_float(char   *line,
		int     n, 
		float *return_value){
  int  error=ERROR_NONE;
  char temp_char[2],temp_char_old[2];
  int  j,k,flag=FALSE;
  strcpy(temp_char_old," ");
  for(k=0,j=0;j<strlen(line);j++) {
    strncpy(temp_char,&(line[j]),sizeof(char));
    sprintf(temp_char,"%c",line[j]);
    if(strcmp(temp_char," ")) {
      if(!strcmp(temp_char_old," ")) {
	k++;
	if(k==n){
	  sscanf(&(line[j]),"%f",return_value);
	  flag=TRUE;
	  j=strlen(line);
	}
      }
    }
    strcpy(temp_char_old,temp_char);
  }
  if(!flag)
    error=ERROR_LINE_TOO_SHORT;
  return(error);
}
