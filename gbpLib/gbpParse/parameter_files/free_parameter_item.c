#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_parameter_files.h>

void free_parameter_item(parameter_item_info **param_item){
   SID_free(SID_FARG (*param_item)->data);
   SID_free(SID_FARG (*param_item));
}

