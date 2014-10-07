#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

int fetch_parameter_list_item(parameter_list_info *param_list,const char *name,parameter_item_info **item){
   (*item)=NULL;
   parameter_item_info *current=param_list->first;
   while(current!=NULL && (*item)==NULL){
      if(!strcmp(current->name,name))
         (*item)=current;
      current=current->next;
   }
   return((*item)!=NULL);
}

