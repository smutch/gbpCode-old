#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMisc.h>

void init_array_log(double   val_min,
                    double   val_max,
                    int      n_val,
                    double **val){
  int    i;
  double step;
  if(n_val>1){
    step=pow(val_max/val_min,1./((double)(n_val-1)));
    if((*val)==NULL)
       (*val)=(double *)SID_malloc(sizeof(double)*n_val);
    (*val)[0]=val_min;
    for(i=1;i<n_val-1;i++)
      (*val)[i]=(*val)[i-1]*step;
    (*val)[n_val-1]=val_max;
  }
  else
    fprintf(stderr,"n<2 in init_array_linear. Array not set.\n");
}

