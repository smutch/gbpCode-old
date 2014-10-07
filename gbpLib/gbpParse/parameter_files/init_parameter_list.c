#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

void init_parameter_list(parameter_list_info **param_list){
   if((*param_list)==NULL)
      (*param_list)=(parameter_list_info *)SID_malloc(sizeof(parameter_list_info));
   (*param_list)->first=NULL;
   (*param_list)->last =NULL;
   (*param_list)->count=0;
}

