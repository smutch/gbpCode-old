#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>

void strip_path(char *string){
  int  i_char;
  int  j_char;
  int  i_start;
  int  string_size;
  char temp_char[2];

  string_size=strlen(string);
  for(i_char=0,i_start=0;i_char<string_size;i_char++){
    strncpy(temp_char,&(string[i_char]),1);
    sprintf(temp_char,"%c",string[i_char]);
    if(!strcmp(temp_char,"/")) 
      i_start=i_char+1;
  }
  if(i_start>0){
    for(i_char=0,j_char=i_start;j_char<string_size;i_char++,j_char++)
      strncpy(&(string[i_char]),&(string[j_char]),1);
    sprintf(&(string[i_char]),"\0");
  }
}

