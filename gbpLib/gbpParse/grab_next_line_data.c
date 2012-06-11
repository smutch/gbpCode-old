#if USE_GETLINE==0
  #if USE_MPI==0
    #define  _GNU_SOURCE
  #endif
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse.h>
int grab_next_line_data(FILE *fp,char **line, size_t *n){
  int rval;
  rval=getline(line,n,fp);
  while(!feof(fp) && check_comment(*line))
    rval=getline(line,n,fp);
  return(rval);
}
