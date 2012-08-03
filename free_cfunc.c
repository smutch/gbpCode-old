#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void free_cfunc(cfunc_info *cfunc){
  SID_log("Freeing correllation function...",SID_LOG_OPEN);
  free_cosmo(&(cfunc->cosmo));

  int i_jack;
  for(i_jack=0;i_jack<=cfunc->n_jack_total;i_jack++){
    SID_free(SID_FARG cfunc->DD_l1D[i_jack]);
    SID_free(SID_FARG cfunc->DR_l1D[i_jack]);
    SID_free(SID_FARG cfunc->RR_l1D[i_jack]);
    SID_free(SID_FARG cfunc->DD_1D[i_jack]);
    SID_free(SID_FARG cfunc->DR_1D[i_jack]);
    SID_free(SID_FARG cfunc->RR_1D[i_jack]);
    SID_free(SID_FARG cfunc->DD_2D[i_jack]);
    SID_free(SID_FARG cfunc->DR_2D[i_jack]);
    SID_free(SID_FARG cfunc->RR_2D[i_jack]);
  }
  SID_free(SID_FARG cfunc->DD_l1D);
  SID_free(SID_FARG cfunc->DR_l1D);
  SID_free(SID_FARG cfunc->RR_l1D);
  SID_free(SID_FARG cfunc->DD_1D);
  SID_free(SID_FARG cfunc->DR_1D);
  SID_free(SID_FARG cfunc->RR_1D);
  SID_free(SID_FARG cfunc->DD_2D);
  SID_free(SID_FARG cfunc->DR_2D);
  SID_free(SID_FARG cfunc->RR_2D);

  int i_run;
  for(i_run=0;i_run<4;i_run++){
     SID_free(SID_FARG cfunc->CFUNC_l1D[i_run]);
     SID_free(SID_FARG cfunc->dCFUNC_l1D[i_run]);
     SID_free(SID_FARG cfunc->COVMTX_l1D[i_run]);
     SID_free(SID_FARG cfunc->CFUNC_1D[i_run]);
     SID_free(SID_FARG cfunc->dCFUNC_1D[i_run]);
     SID_free(SID_FARG cfunc->COVMTX_1D[i_run]);
     SID_free(SID_FARG cfunc->CFUNC_2D[i_run]);
     SID_free(SID_FARG cfunc->dCFUNC_2D[i_run]);
     SID_free(SID_FARG cfunc->COVMTX_2D[i_run]);
  }
  SID_free(SID_FARG cfunc->CFUNC_l1D);
  SID_free(SID_FARG cfunc->dCFUNC_l1D);
  SID_free(SID_FARG cfunc->COVMTX_l1D);
  SID_free(SID_FARG cfunc->CFUNC_1D);
  SID_free(SID_FARG cfunc->dCFUNC_1D);
  SID_free(SID_FARG cfunc->COVMTX_1D);
  SID_free(SID_FARG cfunc->CFUNC_2D);
  SID_free(SID_FARG cfunc->dCFUNC_2D);
  SID_free(SID_FARG cfunc->COVMTX_2D);

  SID_log("Done.",SID_LOG_CLOSE);
}

