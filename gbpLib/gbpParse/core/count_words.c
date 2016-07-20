#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse_core.h>

int count_words(char   *line){
  int  error=ERROR_NONE;
  char temp_char[2],temp_char_old[2];
  int  j,flag=FALSE;
  int  n_words;
  strcpy(temp_char_old," ");
  for(n_words=0,j=0;j<strlen(line);j++) {
    strncpy(temp_char,&(line[j]),1);
    sprintf(temp_char,"%c",line[j]);
    if(strcmp(temp_char," ")) {
      if(!strcmp(temp_char_old," "))
	n_words++;
    }
    strcpy(temp_char_old,temp_char);
  }
  return(n_words);
}
