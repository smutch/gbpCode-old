#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

void prep_profile_custom(double  *ordinate,
			int      n_data,
			double  *ordinate_min,
			double  *ordinate_max,
			int      n_bin,
			int    **bin_id){
  int     i,j,k;
  
  fprintf(stderr,"  Preparing profile bins...");

  (*bin_id)=(int *)malloc(sizeof(int)*n_data);
  for(i=0;i<n_data;i++)
    (*bin_id)[i]=-1;

  for(i=0;i<n_data;i++){
    j=0;
    while((*bin_id)[i]<0 && j<n_bin){
      if(ordinate[i]>=ordinate_min[j] &&
	 ordinate[i]< ordinate_max[j])
	(*bin_id)[i]=j;
      else
	j++;
    }
  }
}
