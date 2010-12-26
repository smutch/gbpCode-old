#include <math.h>
#include <gbpLib.h>
#include <gbpFFT.h>

double k_mag_field_FFT(field_info *FFT,int *i_k){
  int    i_d,j_d;
  double k_mag=0.;
  for(j_d=0;j_d<FFT->n_d;j_d++){
    i_d=j_d;
#ifdef USE_MPI
    if(j_d==1)      i_d=0;
    else if(j_d==0) i_d=1;
#endif
    k_mag+=FFT->k_field[i_d][i_k[i_d]]*
           FFT->k_field[i_d][i_k[i_d]];
  }
  return(sqrt(k_mag));
}

