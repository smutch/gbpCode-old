#include <gbpLib.h>
#include <gbpMisc.h>
void force_periodic(GBPREAL *coord,GBPREAL min,GBPREAL box_size){
  if((*coord)<min)
    (*coord)+=box_size;
  else if((*coord)>=(min+box_size))
    (*coord)-=box_size;
}
