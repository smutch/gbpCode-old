#if USE_GETLINE==0
  #if USE_MPI==0
    #ifndef _GNU_SOURCE
       #define  _GNU_SOURCE
    #endif
  #endif
#endif
#include <stdio.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>

int count_lines_parameters(FILE *fp){
  int     n_lines=0;
  char   *line   =NULL;
  size_t  n      =0;
  int     r;
  while(!feof(fp)){
    r=getline(&line,&n,fp);
    if(r>0 && !check_parameter(line)) 
      n_lines++;
  }
  SID_free(SID_FARG line);
  rewind(fp);
  return(n_lines);
}
