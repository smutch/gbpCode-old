#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

void make_histogram(double  *data,
		    int      n_data,
		    int      n_bin,
		    double  *bin_min,
		    double  *bin_max,
		    int    **histogram){
  int     i,j;

  fprintf(stderr,"  Making histogram...");

  (*histogram)=(int *)malloc(sizeof(int)*n_bin);
  for(i=0;i<n_bin;i++)
    (*histogram)[i]=0;
  for(i=0;i<n_data;i++){
    for(j=0;j<n_bin;j++){
      if(data[i]>=bin_min[j] && data[i]<bin_max[j]){
	(*histogram)[j]++;
	j=n_bin;
      }
    }
  }
  fprintf(stderr,"Done.\n");
}
