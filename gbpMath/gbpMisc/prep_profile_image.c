#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <common.h>

void prep_profile_image(double  *image,
 	                int      n_x,
	                int      n_y,
	                double   x_cen,
	                double   y_cen,
	                int      n_bin,
	                double   first_bin_min_size,
	                double   min_amount,
	                double   r_max_max,
	                int      bin_type,
	                int    **bin_image,
	                double **r_min,
	                double **r_max,
	                int     *n_bin_out){
  int     n_pix;
  double *r_pix;
  int    *order;
  int     i,j,k;
  double  R_next;
  double  R_step;
  double  r_pix_next;
  double  accumulator;
  double  R_max;
  
  fprintf(stderr,"  Creating profile bins...");
  n_pix=n_x*n_y;
  r_pix=(double *)malloc(sizeof(double)*n_pix);
  for(j=0,k=0;j<n_y;j++){
    for(i=0;i<n_x;i++){
      r_pix[k]=pow(((double)i-x_cen)*((double)i-x_cen)+
		   ((double)j-y_cen)*((double)j-y_cen),0.5);
      k++;
    }
  }
  sort_d(r_pix,n_pix,&order,TRUE);
  R_max=MIN(r_max_max,r_pix[n_pix-1]);

  if(bin_type==LOG_BIN){
    R_next=MAX(first_bin_min_size,r_pix[MIN(1,n_pix-1)]);
    R_step=
      pow(R_max/R_next,1.0/(double)(n_bin-1));
  }
  else{
    R_step=
      (R_max/(double)(n_bin));    
    R_next=R_step;
    if(R_next<first_bin_min_size ||
       R_step<r_pix[MIN(1,n_pix-1)]){
      R_next=MAX(first_bin_min_size,r_pix[MIN(1,n_pix-1)]);
      R_step=
	((R_max-R_next)/(double)(n_bin-1));    
    }
  }

  (*bin_image)=(int *)malloc(sizeof(int)*n_pix);
  for(i=0;i<n_pix;i++)
    (*bin_image)[i]=-1;
  (*r_min)    =(double *)malloc(sizeof(double)*n_bin);
  (*r_max)    =(double *)malloc(sizeof(double)*n_bin);
  for(i=0,j=0,(*r_min)[0]=0.0,accumulator=0.0;
      i<n_pix && j<n_bin && r_pix[i]<=R_max;
      i++){
    r_pix_next=r_pix[MIN(i+1,n_pix-1)];
    if(r_pix[i]>=R_next && 
       r_pix_next>r_pix[i]){
      if(accumulator>min_amount){
	if(bin_type==LOG_BIN) R_next*=R_step;
	else                  R_next+=R_step;
	j++;
	if(j<n_bin)
	  (*r_min)[j]=(*r_max)[j-1];
	accumulator=0.0;
      }
      else {
	R_next=r_pix_next;
	if(bin_type==LOG_BIN){
	  R_step=
	    pow(R_max/R_next,1.0/(double)MAX(1,n_bin-j-1));
	}
	else if(bin_type==LINEAR_BIN){
	  R_step=
	    ((R_max-R_next)/(double)(n_bin-j));    
	}
      }
    }
    (*bin_image)[order[i]]=j;
    (*r_max)[j] =r_pix[i];
    accumulator+=image[order[i]];
  }
  *n_bin_out=MIN(n_bin,j+1);
  fprintf(stderr,"%d bins created...Done.\n",(*n_bin_out));
  free(r_pix);
  free(order);
}
