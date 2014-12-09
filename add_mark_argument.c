#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void add_mark_argument(render_info *render,const char *species,int value,const char *type,...){
   va_list vargs;
   va_start(vargs,type);
   mark_arg_info *new_arg;
   create_mark_argument(render,&new_arg);
   strcpy(new_arg->species,species);
   new_arg->value=(char)value;
   strcpy(new_arg->type,type);
   if(!strcmp(new_arg->type,"sphere")){
      for(int i_val=0;i_val<4;i_val++)
         new_arg->dval[i_val]=va_arg(vargs,double);
   }
   else if(!strcmp(new_arg->type,"group_index")){
      for(int i_val=0;i_val<1;i_val++)
         new_arg->ival[i_val]=va_arg(vargs,int);
   }
   else if(!strcmp(new_arg->type,"subgroup_index")){
      for(int i_val=0;i_val<1;i_val++)
         new_arg->ival[i_val]=va_arg(vargs,int);
   }
   else if(!strcmp(new_arg->type,"tree_id")){
      for(int i_val=0;i_val<1;i_val++)
         new_arg->ival[i_val]=va_arg(vargs,int);
   }
   else if(!strcmp(new_arg->type,"") || !strcmp(new_arg->type,"*"))
      sprintf(new_arg->type,"all");
   else 
      SID_trap_error("Invalid mark type {%s} in add_mark_argument().",ERROR_LOGIC,new_arg->type);
   if(render->mark_arg_first==NULL)
      render->mark_arg_first=new_arg;
   else
      render->mark_arg_last->next=new_arg;
   render->mark_arg_last=new_arg;
   va_end(vargs);
}

