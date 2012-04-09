#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMisc.h>

int count_words(char   *line){
  int  error=ERROR_NONE;
  char temp_char[2],temp_char_old[2];
  int  j,k,flag=FALSE;
  int  n;
  strcpy(temp_char_old," ");
  for(n=0,j=0;j<strlen(line);j++) {
    strncpy(temp_char,&(line[j]),sizeof(char));
    sprintf(temp_char,"%c\0",line[j]);
    if(strcmp(temp_char," ")) {
      if(!strcmp(temp_char_old," ")) {
	n++;
      }
    }
    strcpy(temp_char_old,temp_char);
  }
  return(n);
}
