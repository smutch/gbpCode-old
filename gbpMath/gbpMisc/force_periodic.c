#include <gbpLib.h>
void force_periodic(REAL *coord,REAL min,REAL box_size){
  if((*coord)<min)
    (*coord)+=box_size;
  else if((*coord)>=(min+box_size))
    (*coord)-=box_size;
}
