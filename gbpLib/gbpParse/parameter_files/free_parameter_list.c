#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

int free_parameter_list(parameter_list_info **param_list){
   parameter_item_info *current=(*param_list)->first;
   while(current!=NULL){
      parameter_item_info *next=current->next;
      free_parameter_item(&current);
      current=next;
   }
   SID_free(SID_FARG (*param_list));
}

