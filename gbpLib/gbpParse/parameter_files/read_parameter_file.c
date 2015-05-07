#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>
#include <gbpParse_parameter_files.h>

void read_parameter_file(const char          *filename_in,
                         parameter_list_info *param_list){
   FILE   *fp_in=NULL;
   char   *line =NULL;
   size_t  line_length=0;
   if((fp_in=fopen(filename_in,"r"))==NULL)
      SID_trap_error("Could not open parameter file {%s}.",ERROR_LOGIC,filename_in);
   while(grab_next_line_parameter(fp_in,&line,&line_length)){
      char                 parameter_name[256];
      double               parameter_value;
      parameter_item_info *parameter_item;
      remove_parameter_character(line);
      grab_word(line,1,parameter_name);
      if(!fetch_parameter_list_item(param_list,parameter_name,&parameter_item))
         SID_trap_error("Paramter {%s} is not in the file specification for {%s}.",ERROR_LOGIC,parameter_name,filename_in);
      if(parameter_item->data_type!=SID_CHAR && count_words(line)!=2)
         SID_trap_error("Invalid number of words in parameter given by {%s}.",ERROR_LOGIC,line);
      if(parameter_item->data_type==SID_DOUBLE)
         grab_double(line,2,(double *)(parameter_item->data));
      else if(parameter_item->data_type==SID_INT)
         grab_int(line,2,(int *)(parameter_item->data));
      else if(parameter_item->data_type==SID_SIZE_T)
         grab_size_t(line,2,(size_t *)(parameter_item->data));
      else if(parameter_item->data_type==SID_CHAR)
         grab_tail(line,2,(char *)(parameter_item->data));
      else
         SID_trap_error("Unsupported parameter datatype specified for parameter {%s}.",ERROR_LOGIC,parameter_name);
      parameter_item->flag_set=TRUE;
      parameter_item->n_read++;
   }
   SID_free(SID_FARG line);
   fclose(fp_in);

   // Check that all manditory parameters were read
   parameter_item_info *current=param_list->first;
   while(current!=NULL){
      if(check_mode_for_flag(current->mode,PARAMETER_MODE_MANDITORY) && current->n_read<=0)
         SID_trap_error("Manditory parameter {%s} was not specified in {%s}.",ERROR_LOGIC,current->name,filename_in);
      if(check_mode_for_flag(current->mode,PARAMETER_MODE_UNIQUE) && current->n_read>1)
         SID_trap_error("Parameter {%s} was specified %d times in {%s}.",ERROR_LOGIC,current->name,current->n_read,filename_in);
      current=current->next;
   }
}

