#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

void make_profile(double  *data,
                  double  *weight,
                  int     *bin_id,
		  int      n_data,
		  int      n_bin,
                  int      avg_flag,
		  double **profile){
  int     i,bin;
  int    *n_pix_bin;
  double *norm;

  fprintf(stderr,"  Making profile...");

  (*profile)=(double *)malloc(sizeof(double)*n_bin);
  norm      =(double *)malloc(sizeof(double)*n_bin);
  for(i=0;i<n_bin;i++){
    (*profile)[i]=0.0;
    norm[i]      =0.0;
  }
  for(i=0;i<n_data;i++){
    if(bin_id[i]>=0){
      if(avg_flag!=AVG_WEIGHT){
        (*profile)[bin_id[i]]+=data[i];
        norm[bin_id[i]]+=1.0;
      }
      else{
        (*profile)[bin_id[i]]+=data[i]*weight[i];
        norm[bin_id[i]]+=weight[i];
      }
    }
  }
  if(avg_flag!=AVG_SUM)
    for(i=0;i<n_bin;i++)
      (*profile)[i]/=norm[i];
  free(norm);
  fprintf(stderr,"Done.\n");
}
