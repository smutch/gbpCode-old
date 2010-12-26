#include <common.h>

double calc_median_marked(double *data,
                          int    *marked,
		          int     n_data){
  int     n_marked;
  int    *sort;
  int     i_mid_1;
  int     i_mid_2;
  int     index_1;
  int     index_2;
  int     i;
  double  rval;
  for(n_marked=0,i=0;i<n_data;i++)
    if(marked[i])
      n_marked++;
  if(n_marked>0){
    sort_d(data,n_data,&sort,FALSE);
    i_mid_1=1+n_marked/2;
    if(n_marked%2)
      i_mid_2=i_mid_1;
    else
      i_mid_2=i_mid_1+1;
    i=0;
    index_1=0;
    while(i<i_mid_1){
      if(marked[sort[index_1]]) i++;
      if(i<i_mid_1)
        index_1++;
    }
    index_2=index_1;
    if(i_mid_1!=i_mid_2){
      while(i<i_mid_2){
        if(marked[sort[index_2]]) i++;
        if(i<i_mid_2)
          index_2++;
      }
    }
    if(i_mid_1==i_mid_2)
      rval=data[sort[index_1]];
    else
      rval=ONE_HALF*(data[sort[index_1]]+data[sort[index_2]]);
    free(sort);
  }
  else
    rval=0.;
  return(rval);
}
