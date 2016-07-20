#if USE_GETLINE==0
  #if USE_MPI==0
    #ifndef _GNU_SOURCE
      #define  _GNU_SOURCE
    #endif
  #endif
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>
int grab_next_line_parameter(FILE *fp,char **line, size_t *n){
  int rval=FALSE;
  if(!feof(fp)){
     rval=getline(line,n,fp);
     while(!feof(fp) && !check_parameter(*line))
       rval=getline(line,n,fp);
     if(!check_parameter(*line))
        rval=FALSE;
  }
  return(rval);
}
