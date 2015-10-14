#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpTrees_build.h>

int tree_case_flags_text(int match_type,const char *separator_string,char **return_string){
  // Allocate the string if necessary
  if((*return_string)==NULL){
     int max_size      =0;
     int separator_size=strlen(separator_string);
     for(int i_parse=0;i_parse<n_tree_case_flag_list;i_parse++)
        max_size+=(strlen(tree_case_flag_list_text[i_parse])+separator_size);
     (*return_string)=(char *)SID_malloc(sizeof(char)*(max_size+2));
  }
  int count=0;
  sprintf((*return_string),"");
  for(int i_parse=0;i_parse<n_tree_case_flag_list;i_parse++){
     if(check_mode_for_flag(match_type,tree_case_flag_list[i_parse])){
        count++;
        if(count==1)
           sprintf((*return_string),"%s",tree_case_flag_list_text[i_parse]);
        else
           sprintf((*return_string),"%s%s%s",(*return_string),separator_string,tree_case_flag_list_text[i_parse]);
     }
  }
  return(count);
}

