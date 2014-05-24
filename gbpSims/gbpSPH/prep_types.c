#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void prep_types(char ***pname, int n_type, ... ){
  va_list  ap;
  int      i;
  int      index;
  char    *name;
  va_start(ap,n_type);
  (*pname)=(char **)SID_malloc(sizeof(char *)*N_GADGET_TYPE);
  for(i=0;i<N_GADGET_TYPE;i++){
    (*pname)[i]=(char *)SID_malloc(sizeof(char)*256);
    strcpy((*pname)[i],"unused");
  }
  for(i=0;i<n_type;i++){
    index=(int)va_arg(ap,int);
    name =(char *)va_arg(ap,char *);
    strcpy((*pname)[index],name);
  }
  va_end(ap);
}
