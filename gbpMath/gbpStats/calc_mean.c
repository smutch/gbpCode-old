#include <gbpLib.h>

double calc_mean(void  *data,
		 size_t n_data,
                 SID_Datatype type){
  double  mean;
  if(n_data<1)
    mean=0.;
  else
    mean=calc_sum(data,n_data,type)/(double)n_data;
  return(mean);
}
