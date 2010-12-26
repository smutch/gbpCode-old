#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
size_t find_index(size_t *y,size_t y_find,size_t  n,size_t *index){
  size_t x_lo;
  size_t x_hi;
  size_t x_mid;
  size_t y_lo;
  size_t y_hi;
  size_t y_mid;
  size_t y_max;
  x_lo=0;
  x_hi=n-1;
  if(index==NULL){
    y_lo=y[x_lo];
    y_hi=y[x_hi];
  }
  else{
    y_lo=y[index[x_lo]];
    y_hi=y[index[x_hi]];
  }
  y_max=y_hi;
  do{
    x_mid=(x_lo+x_hi)/2;
    if(index==NULL)
      y_mid=y[x_mid];
    else
      y_mid=y[index[x_mid]];
    if(y_find<y_mid)
      x_hi=x_mid;
    else if(y_find>y_mid)
      x_lo=x_mid;
    else{
      x_lo=x_mid;
      x_hi=x_mid;
    }
  } while(x_hi-x_lo>2);

  // Make sure we are at the beginning of the wanted
  //  index or, if it does not exist, return the
  //  position of the previous index
  if(x_mid<n-1){
    x_hi=x_mid+1;
    if(index==NULL)
      y_hi=y[x_hi];
    else
      y_hi=y[index[x_hi]];
  }
  while(y_find<=y_mid && x_mid>0){
    y_hi=y_mid;
    x_hi=x_mid;
    x_mid--;
    if(index==NULL)
      y_mid=y[x_mid];
    else
      y_mid=y[index[x_mid]];
  }
  if(y_find==y_mid) // needed for first element
    return(x_mid);
  else if(y_find==y_hi)
    return(x_hi);
  else if(y_find==y_max)
    return(n-1);
  else
    return(x_mid);
}
