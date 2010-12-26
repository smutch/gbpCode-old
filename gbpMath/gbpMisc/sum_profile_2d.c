#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

void sum_profile_2d(double  *r_max,
		    double  *data,
		    int      n_bin,
		    double **profile){
  int     i;
  double  area;

  (*profile)=(double *)malloc(sizeof(double)*n_bin);
  for(i=0;i<n_bin;i++){
    if(i==0){
      area=PI*r_max[i]*r_max[i];
      (*profile)[i]=data[i]*area;
    }
    else{
      area=TWO_PI*r_max[i]*(r_max[i]-r_max[i-1]);
      (*profile)[i]=(*profile)[i-1]+data[i]*area;
    }
  }

}
