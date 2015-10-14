#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void seal_render(render_info *render){
  scene_info *current;
  scene_info *next;
  SID_log("Sealing render...",SID_LOG_OPEN);
  seal_scenes(render->scenes);
  seal_render_camera(render);
  render->n_frames=render->last_scene->last_frame+1;
  render->sealed=TRUE;

  SID_log("Done.",SID_LOG_CLOSE);
}

