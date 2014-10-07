#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

int fetch_parameter_data(parameter_list_info *param_list,const char *name,void *value){
   parameter_item_info *item;
   if(!fetch_parameter_list_item(param_list,name,&item))
      return(FALSE);
   else if(!(item->flag_set))
      return(FALSE);
   else{
      memcpy(value,item->data,item->data_size);
      return(TRUE);
   }
}

