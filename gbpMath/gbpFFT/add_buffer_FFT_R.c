#include <gbpLib.h>
#include <gbpFFT.h>

void add_buffer_FFT_R(field_info *FFT){
  size_t index_pad;
  size_t i_FFT;
  if(!(FFT->flag_padded)){
    SID_log("Adding FFTW padding...",SID_LOG_OPEN);
    for(i_FFT=FFT->n_field_R_local-1;i_FFT>0;i_FFT--){
      index_pad=pad_index_FFT_R(FFT,i_FFT);
//if(index_pad>=FFT->total_local_size) fprintf(stderr,"Arg! %d %lld %lld %lld\n",SID.My_rank,i_FFT,index_pad,FFT->total_local_size);
      if(index_pad!=i_FFT)
        FFT->field_local[index_pad]=FFT->field_local[i_FFT];
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  FFT->flag_padded=TRUE;
}

