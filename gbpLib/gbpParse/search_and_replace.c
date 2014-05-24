#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse.h>

// TODO: Need to deal with replace_length!=search_length properly.

int search_and_replace(char *string,const char *search,const char *replace){
  int string_length;
  int search_length;
  int replace_length;
  int i_string;
  int n_replace;

  string_length =strlen(string);
  search_length =strlen(search);
  replace_length=strlen(replace);

  if(search_length!=replace_length) 
     SID_log_warning("search_length!=replace_length in search_and_replace()",SID_WARNING_DEFAULT);

  for(i_string=0,n_replace=0;i_string<string_length-search_length; i_string++){
    if(!strncmp(&(string[i_string]),search,search_length)){
      strncpy(&(string[i_string]),replace,replace_length);
      n_replace++;
    }
  }
  return(n_replace);
}

