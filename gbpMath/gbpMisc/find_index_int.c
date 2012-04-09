#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpMisc.h>
int find_index_int(int *y,int y_find,int  n,size_t *index){
  int x_lo;
  int x_hi;
  int x_mid;
  int y_lo;
  int y_hi;
  int y_mid;
  int y_max;
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
