#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void create_mark_argument(render_info *render,mark_arg_info **new_arg){
   (*new_arg)=(mark_arg_info *)SID_malloc(sizeof(mark_arg_info));
   sprintf((*new_arg)->species,"");
   sprintf((*new_arg)->type,   "");
   (*new_arg)->next=NULL;
}

