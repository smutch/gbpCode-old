#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void set_render_scale(render_info *render,double RGB_min,double RGB_max,double Y_min,double Y_max,double Z_min,double Z_max){
  render->camera->RGB_range[0]=RGB_min;
  render->camera->RGB_range[1]=RGB_max;
  render->camera->Y_range[0]  =Y_min;
  render->camera->Y_range[1]  =Y_max;
  render->camera->Z_range[0]  =Z_min;
  render->camera->Z_range[1]  =Z_max;
}
