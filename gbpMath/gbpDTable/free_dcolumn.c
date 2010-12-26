#include <stdio.h>
#include <common.h>
int free_dcolumn(dcolumn **remove){
  free((*remove)->data);
  free((*remove));
  (*remove)=NULL;
  return(ERROR_NONE);
}
