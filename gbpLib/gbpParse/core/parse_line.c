#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>

int parse_line(char *line,
	       int   n_return, ...){ 
  char     temp_char[2],temp_char_old[2];
  int      i_line;
  int      i_word;
  int      i_parse;
  int      flag_done=FALSE;
  void   **return_values;
  int     *return_types;
  int     *return_words;
  va_list  vargs;
  int      error=ERROR_NONE;

  // Build parse list
  va_start(vargs,n_return);
  return_words =(int   *)SID_malloc(sizeof(int   )*n_return);
  return_values=(void **)SID_malloc(sizeof(void *)*n_return);
  return_types =(int   *)SID_malloc(sizeof(int   )*n_return);
  for(i_parse=0;i_parse<n_return;i_parse++){
    return_words[i_parse] =(int)   va_arg(vargs,int);
    return_values[i_parse]=(void *)va_arg(vargs,void *);
    return_types[i_parse] =(int)   va_arg(vargs,int);
  }

  // Perform parsing
  if(n_return>0){
    strcpy(temp_char_old," ");
    for(i_word=0,i_line=0,i_parse=0;i_line<strlen(line);i_line++) {
      strncpy(temp_char,&(line[i_line]),1);
      sprintf(temp_char,"%c\0",line[i_line]);
      if(strcmp(temp_char," ")) {
        // If we have whitespace, then we have a new word
        if(!strcmp(temp_char_old," ")) {
	  i_word++;
	  if(i_word==return_words[i_parse]){
            switch(return_types[i_parse]){
              case SID_DOUBLE:
  	        sscanf(&(line[i_line]),"%lf",(double *)return_values[i_parse]);
                break;
              case SID_LONG:
  	        sscanf(&(line[i_line]),"%ld",(long *)return_values[i_parse]);
                break;
              case SID_FLOAT:
  	        sscanf(&(line[i_line]),"%f",(float *)return_values[i_parse]);
                break;
              case SID_INT:
  	        sscanf(&(line[i_line]),"%d",(int *)return_values[i_parse]);
                break;
              case SID_SIZE_T:
  	        sscanf(&(line[i_line]),"%zu",(size_t *)return_values[i_parse]);
                break;
              default:
                SID_trap_error("Type %d is not supported in parse_line!",ERROR_LOGIC,return_types[i_parse]);
                break;
            }
            i_word++;
            i_parse++;
            if(i_parse==n_return){
	      flag_done=TRUE;
	      break;
            }
 	  }
        }
      }
      strcpy(temp_char_old,temp_char);
    }
    if(!flag_done)
      error=ERROR_LINE_TOO_SHORT;
  }

  // Clean-up
  SID_free(SID_FARG return_words);
  SID_free(SID_FARG return_values);
  SID_free(SID_FARG return_types);
  va_end(vargs);

  return(error);
}
