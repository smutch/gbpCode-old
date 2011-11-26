// Free memory allocated by initialize_field()
#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpDomain.h>

void free_field(field_info *FFT){
  int i_d;

  free(FFT->n_R_local);
  free(FFT->n_k_local);
  free(FFT->i_R_start_local);
  free(FFT->i_k_start_local);
  free(FFT->i_R_stop_local);
  free(FFT->i_k_stop_local);

  // Free FFTs
  #if USE_MPI
    rfftwnd_mpi_destroy_plan(FFT->plan);
    rfftwnd_mpi_destroy_plan(FFT->iplan);
  #else
    rfftwnd_destroy_plan(FFT->plan);
    rfftwnd_destroy_plan(FFT->iplan);
  #endif

  // Free field arrays
  for(i_d=0;i_d<FFT->n_d;i_d++){
    free(FFT->k_field[i_d]);
    free(FFT->R_field[i_d]);
  }
  free(FFT->k_field);
  free(FFT->R_field);
  free(FFT->field_local);
  free(FFT->dR);
  free(FFT->dk);
  free(FFT->k_Nyquist);
  free(FFT->n);
  free(FFT->L);
}

