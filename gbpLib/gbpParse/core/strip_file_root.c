#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse_core.h>

void strip_file_root(char *string){
  int  i_char;
  int  j_char;
  int  flag_continue;
  int  string_start;
  char temp_char[2];

  string_start=strlen(string)-1;
  for(i_char=string_start,flag_continue=TRUE;i_char>=0 && flag_continue;i_char--){
    strncpy(temp_char,&(string[i_char]),1);
    sprintf(temp_char,"%c",string[i_char]);
    if(!strcmp(temp_char,"/")){
      if(i_char<string_start)
         sprintf(&(string[i_char+1]),"\0");
      flag_continue=FALSE;
    }
  }
  if(flag_continue)
    sprintf(string,"./");
}

