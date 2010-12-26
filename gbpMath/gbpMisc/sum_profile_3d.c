#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

void sum_profile_3d(double  *r_max,
                    double  *data,
		    int      n_bin,
		    double **profile){
  int     i;
  double  volume;

  (*profile)=(double *)malloc(sizeof(double)*n_bin);
  for(i=0;i<n_bin;i++){
    if(i==0){
      volume=FOUR_THIRDS_PI*r_max[i]*r_max[i]*r_max[i];
      (*profile)[i]=data[i]*volume;
    }
    else{
      volume=FOUR_PI*r_max[i]*r_max[i]*(r_max[i]-r_max[i-1]);
      (*profile)[i]=(*profile)[i-1]+data[i]*volume;
    }
  }

}
