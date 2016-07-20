#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void init_scene(scene_info **scene){
  SID_log("Initializing scene...",SID_LOG_OPEN);
  (*scene)=(scene_info *)SID_malloc(sizeof(scene_info));
  (*scene)->n_frames         =1;
  (*scene)->n_perspectives   =0;
  (*scene)->perspectives     =NULL;
  (*scene)->first_perspective=NULL;
  (*scene)->last_perspective =NULL;
  (*scene)->evolve           =NULL;
  add_scene_perspective((*scene));
  init_perspective(&((*scene)->evolve),RENDER_INIT_EVOLVE);
  init_perspective_interp(&((*scene)->interp));
  (*scene)->sealed     =FALSE;
  (*scene)->next       =NULL;
  SID_log("Done.",SID_LOG_CLOSE);
}

