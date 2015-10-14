#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void free_mark_arguments(mark_arg_info **argument){
   while((*argument)!=NULL){
      mark_arg_info *next=(*argument)->next;
      SID_free(SID_FARG (*argument));
      (*argument)=next;
   }
}

