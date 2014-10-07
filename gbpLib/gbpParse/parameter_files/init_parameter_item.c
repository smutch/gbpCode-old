#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

void init_parameter_item(parameter_item_info **param_item,
                         const char           *name,
                         SID_Datatype          data_type,
                         int                   mode){
   if((*param_item)==NULL)
      (*param_item)=(parameter_item_info *)SID_malloc(sizeof(parameter_item_info));
   strcpy((*param_item)->name,name);
   (*param_item)->data_type=data_type;
   SID_Type_size(data_type,&((*param_item)->data_size));
   if(data_type==SID_CHAR)
      (*param_item)->data_size*=PARAMETER_STRING_LENGTH;
   (*param_item)->n_read   =0;
   (*param_item)->flag_set =FALSE;
   (*param_item)->mode     =mode;
   (*param_item)->data     =SID_calloc((*param_item)->data_size);
   (*param_item)->next     =NULL;
}

