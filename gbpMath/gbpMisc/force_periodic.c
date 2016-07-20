#include <gbpLib.h>
#include <gbpMisc.h>
#define  N_MAX_LOCAL 10

void force_periodic(GBPREAL *coord,GBPREAL min,GBPREAL box_size){
  int n=0;
  while((*coord)<min){
    if(n>N_MAX_LOCAL) SID_trap_error("N_MAX reached in force_periodic. (1)",ERROR_LOGIC);
    (*coord)+=box_size;
    n++;
  }
  if(n==0){
     while((*coord)>=(min+box_size)){
        if(n>N_MAX_LOCAL) SID_trap_error("N_MAX reached in force_periodic. (2)",ERROR_LOGIC);
        (*coord)-=box_size;
        n++;
     }
  }
}
