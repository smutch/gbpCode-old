#include <stdio.h>
#include <stdlib.h>
#include <common.h>
int find_index_start(int *x,int x_find,int  n){
  int x_lo;
  int x_hi;
  int x_mid;
  int i;

  x_lo =x[0];
  x_hi =x[n-1];
  y_lo =y[0];
  y_hi =y[n-1];
  do{
    x_mid=(x_lo+x_hi)/2;
    if(x_find<x_mid)
      x_hi=x_mid;
    else if(x_find>x_mid)
      x_lo=x_mid;
    else{
      x_lo=x_mid;
      x_hi=x_mid;
    }
  } while(x_hi-x_lo>2);
  return(x_lo);
}
