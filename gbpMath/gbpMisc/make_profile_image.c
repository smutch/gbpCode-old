#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>

void make_profile_image(double  *image,
 		        int     *bin_image,
		        int      n_pix,
		        int      n_bin,
		        double **profile,
		        int      avg_flag){
  int  i,bin;
  int *n_pix_bin;

  (*profile)=(double *)malloc(sizeof(double)*n_bin);
  n_pix_bin =(int    *)malloc(sizeof(int   )*n_bin);
  for(i=0;i<n_bin;i++){
    (*profile)[i]=0.0;
    n_pix_bin[i]     =0;
  }
  for(i=0;i<n_pix;i++){
    bin=bin_image[i];
    if(bin>=0){
      (*profile)[bin]+=image[i];
      n_pix_bin[bin]++;
    }
  }
  if(avg_flag){
    for(i=0;i<n_bin;i++){
      if(n_pix_bin[i]>0)
	(*profile)[i]/=(double)n_pix_bin[i];
    }
  }
  free(n_pix_bin);
}
