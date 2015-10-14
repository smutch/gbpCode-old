#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void copy_perspective(perspective_info *from,perspective_info *to){
  SID_log("Copying perspective...",SID_LOG_OPEN);
  to->p_o[0]       =from->p_o[0];
  to->p_o[1]       =from->p_o[1];
  to->p_o[2]       =from->p_o[2];
  to->radius       =from->radius;
  to->FOV          =from->FOV;
  to->focus_shift_x=from->focus_shift_x;
  to->focus_shift_y=from->focus_shift_y;
  to->theta        =from->theta;
  to->zeta         =from->zeta;
  to->phi          =from->phi;
  to->time         =from->time;
  SID_log("Done.",SID_LOG_CLOSE);
}

