#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse.h>

int grab_word(char  *line,
	      int    n, 
	      char **return_value,
              int   *size){
  int  error=ERROR_NONE;
  char temp_char[2],temp_char_old[2];
  int  j,k,l,size_temp,flag=FALSE;
  strcpy(temp_char_old," ");
  for(j=0,k=0;j<strlen(line) && !flag;j++) {
    sprintf(temp_char,"%c\0",line[j]);
    if(strcmp(temp_char," ")) {
      if(!strcmp(temp_char_old," ")) {
	k++;
	if(k==n){
          l=j;
          while(strcmp(temp_char," ") && j<strlen(line)){
            j++;
            sprintf(temp_char,"%c\0",line[j]);
          }
          size_temp=j-l+1;
          if((*size)<=0)
            *return_value=(char *)malloc(sizeof(char)*size_temp);
          else if((*size)<size_temp){
            *size        =size_temp;
            *return_value=(char *)realloc(*return_value,*size);
          }
          strncpy(*return_value,&(line[l]),(size_temp-1)*sizeof(char));
          sprintf(temp_char,"\0");
          strncpy(&((*return_value)[size_temp-1]),temp_char,sizeof(char));
	  flag=TRUE;
	}
      }
    }
    strcpy(temp_char_old,temp_char);
  }
  if(!flag)
    error=ERROR_LINE_TOO_SHORT;
  return(error);
}
