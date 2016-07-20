#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void init_perspective_interp(perspective_interp_info **perspective_interp){
  SID_log("Initializing perspective interpolation structures...",SID_LOG_OPEN);
  (*perspective_interp)=(perspective_interp_info *)SID_malloc(sizeof(perspective_interp_info));
  (*perspective_interp)->p_o[0]=NULL;
  (*perspective_interp)->p_o[1]=NULL;
  (*perspective_interp)->p_o[2]=NULL;
  (*perspective_interp)->radius=NULL;
  (*perspective_interp)->FOV   =NULL;
  (*perspective_interp)->theta =NULL;
  (*perspective_interp)->zeta  =NULL;
  (*perspective_interp)->phi   =NULL;
  (*perspective_interp)->time  =NULL;
  SID_log("Done.",SID_LOG_CLOSE);
}

