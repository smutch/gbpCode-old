#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse.h>

int grab_char(char *line,
	      int   n, 
	      char *return_value){
  int  error=ERROR_NONE;
  char temp_char[2],temp_char_old[2],temp_string[1000];
  int  j,k,flag=FALSE;
  strcpy(temp_char_old," ");
  for(k=0,j=0;j<strlen(line);j++) {
    strncpy(temp_char,&(line[j]),sizeof(char));
    sprintf(temp_char,"%c\0",line[j]);
    if(strcmp(temp_char," ")) {
      if(!strcmp(temp_char_old," ")) {
	k++;
	if(k==n){
	  sscanf(&(line[j]),"%s ",temp_string);
          sprintf(return_value,"%s\0",temp_string);
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
