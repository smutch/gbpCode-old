#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void init_perspective(perspective_info **perspective,int mode){
  SID_log("Initializing perspective...",SID_LOG_OPEN);
  (*perspective)=(perspective_info *)SID_malloc(sizeof(perspective_info));
  (*perspective)->p_o[0]       = 0.;
  (*perspective)->p_o[1]       = 0.;
  (*perspective)->p_o[2]       = 0.;
  (*perspective)->theta        = 0.;
  (*perspective)->zeta         = 0.;
  (*perspective)->FOV          = 0.;
  (*perspective)->focus_shift_x= 0.;
  (*perspective)->focus_shift_y= 0.;
  if(check_mode_for_flag(mode,RENDER_INIT_EVOLVE)){
    (*perspective)->radius= 0.;
    (*perspective)->phi   = 0.;
    (*perspective)->time  = 0.;
  }
  else{
    (*perspective)->radius= 1.;
    (*perspective)->phi   = 1.;
    (*perspective)->time  = 1.;
  }
  (*perspective)->next  = NULL;
  SID_log("Done.",SID_LOG_CLOSE);
}

