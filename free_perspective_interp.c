#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void free_perspective_interp(perspective_interp_info **perspective_interp){
  SID_log("Freeing perspective interpolation structures...",SID_LOG_OPEN);
  free_interpolate(SID_FARG (*perspective_interp)->p_o[0],NULL);
  free_interpolate(SID_FARG (*perspective_interp)->p_o[1],NULL);
  free_interpolate(SID_FARG (*perspective_interp)->p_o[2],NULL);
  free_interpolate(SID_FARG (*perspective_interp)->theta, NULL);
  free_interpolate(SID_FARG (*perspective_interp)->FOV,   NULL);
  free_interpolate(SID_FARG (*perspective_interp)->radius,NULL);
  free_interpolate(SID_FARG (*perspective_interp)->zeta,  NULL);
  free_interpolate(SID_FARG (*perspective_interp)->phi,   NULL);
  free_interpolate(SID_FARG (*perspective_interp)->time,  NULL);
  SID_free(SID_FARG (*perspective_interp));
  SID_log("Done.",SID_LOG_CLOSE);
}

