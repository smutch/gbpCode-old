#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

int prep_profile(double  *ordinate,
	         int      n_data,
	         int      n_bin,
	         double   first_bin_min_size,
	         int      min_number,
	         double   ordinate_max_max,
	         int      bin_type,
	         int    **bin_id,
	         double **ordinate_min,
	         double **ordinate_max,
	         double **ordinate_avg,
	         int     *n_bin_out){
  int    *ordinate_id;
  int     i,j,k;
  double  next_ordinate;
  double  ordinate_next;
  double  ordinate_step;
  int     accumulator;
  double  ordinate_stop;
  int     status=ERROR_NONE;
  
  fprintf(stderr,"  Preparing profile bins...");

  sort_d(ordinate,n_data,&ordinate_id,FALSE);

  ordinate_stop=MIN(ordinate_max_max,ordinate[ordinate_id[n_data-1]]);

  if(bin_type==LOG_BIN){
    ordinate_next=MAX(first_bin_min_size,ordinate[ordinate_id[MIN(1,n_data-1)]]);
    ordinate_step=
      pow(ordinate_stop/ordinate_next,1.0/(double)(n_bin-1));
  }
  else{
    ordinate_step=
      (ordinate_stop/(double)(n_bin));    
    ordinate_next=ordinate_step;
    if(ordinate_next<first_bin_min_size ||
       ordinate_step<ordinate[ordinate_id[MIN(1,n_data-1)]]){
      ordinate_next=MAX(first_bin_min_size,ordinate[ordinate_id[MIN(1,n_data-1)]]);
      ordinate_step=
	((ordinate_stop-ordinate_next)/(double)(n_bin-1));    
    }
  }

  (*bin_id)=(int *)malloc(sizeof(int)*n_data);
  for(i=0;i<n_data;i++)
    (*bin_id)[i]=-1;
  (*ordinate_min)    =(double *)malloc(sizeof(double)*n_bin);
  (*ordinate_max)    =(double *)malloc(sizeof(double)*n_bin);
  (*ordinate_avg)    =(double *)malloc(sizeof(double)*n_bin);
  for(i=0;i<n_bin;i++)
    (*ordinate_avg)[i]=0.0;
  for(i=0,j=0,(*ordinate_min)[0]=0.0,accumulator=0;
      i<n_data && j<n_bin && ordinate[ordinate_id[i]]<=ordinate_stop;
      i++){
    next_ordinate=ordinate[ordinate_id[MIN(i+1,n_data-1)]];
    if(ordinate[ordinate_id[i]]>=ordinate_next && 
       next_ordinate>ordinate[ordinate_id[i]]){
      if(accumulator>min_number){
	if(bin_type==LOG_BIN) ordinate_next*=ordinate_step;
	else                  ordinate_next+=ordinate_step;
        (*ordinate_avg)[j]/=(double)accumulator;
	j++;
	if(j<n_bin)
	  (*ordinate_min)[j]=(*ordinate_max)[j-1];
	accumulator=0;
      }
      else {
	ordinate_next=next_ordinate;
	if(bin_type==LOG_BIN){
	  ordinate_step=
	    pow(ordinate_stop/ordinate_next,1.0/(double)MAX(1,n_bin-j-1));
	}
	else if(bin_type==LINEAR_BIN){
	  ordinate_step=
	    ((ordinate_stop-ordinate_next)/(double)(n_bin-j));    
	}
      }
    }
    (*ordinate_avg)[j]+=ordinate[ordinate_id[i]];
    (*bin_id)[ordinate_id[i]]=j;
    (*ordinate_max)[j]       =ordinate[ordinate_id[i]];
    accumulator++;
  }
  for(i=0;i<n_data;i++)
    if((*bin_id)[i]==j)
      (*bin_id)[i]=-1;
  *n_bin_out=MIN(n_bin,j);
  fprintf(stderr,"%d bins created...Done.\n",(*n_bin_out));
  free(ordinate_id);
  return(status);
}
