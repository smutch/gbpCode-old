#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpParse_parameter_files.h>

int remove_parameter_character(char *line){
   char temp_char[2];
   int  j,flag=TRUE,rval=FALSE;
   for(j=0;j<strlen(line) && flag;j++){
      sprintf(temp_char,"%c",line[j]);
      if(strcmp(temp_char," ")){
         if(!strcmp(temp_char,GBPPARSE_PARAMETER_CHARACTER)){
            char replace[1];
            sprintf(replace," ");
            strncpy(&(line[j]),replace,1);
            rval=TRUE;
            flag=FALSE;
         }
      }
   }
   return(rval);
}

