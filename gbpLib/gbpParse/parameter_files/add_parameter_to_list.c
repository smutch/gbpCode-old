#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

void add_parameter_to_list(parameter_list_info *param_list,
                           const char          *name,
                           SID_Datatype         data_type,
                           int                  mode){
   parameter_item_info *new_item=NULL;
   init_parameter_item(&new_item,name,data_type,mode);
   if(param_list->first==NULL)
      param_list->first=new_item;
   else
      param_list->last->next=new_item;
   param_list->last=new_item;
   param_list->count++;
}

