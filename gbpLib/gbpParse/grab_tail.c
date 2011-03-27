#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse.h>

int grab_tail(char *line,
	      int   n, 
	      char *return_value){
  int  error=ERROR_NONE;
  char temp_char[2],temp_char_old[2];
  int  j,k,l,m,flag=FALSE;
  strcpy(temp_char_old," ");
  for(k=0,j=0;j<strlen(line);j++) {
    strncpy(temp_char,&(line[j]),1);
    sprintf(temp_char,"%c",line[j]);
    if(strcmp(temp_char," ")) {
      if(!strcmp(temp_char_old," ")) {
	k++;
	if(k==n){
          for(l=-1,m=strlen(line)-1;m>=j;m--){
            sprintf(temp_char,"%c",line[m]);
            if(l<0 && !strcmp(temp_char_old," "))
              l=m;
          }
	  strncpy(return_value,&(line[j]),(l-j));
          strcpy(&(return_value[(l-j)]),"\0");
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
