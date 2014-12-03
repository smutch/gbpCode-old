#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void add_render_scene(render_info *render){
  scene_info *current_s;
  scene_info *last_s;
  scene_info *next_s;
  scene_info *new_s;
  SID_log("Adding render scene...",SID_LOG_OPEN);
  // Seal the previous scene
  seal_scenes(render->scenes);
  // Create a new scene
  init_scene(&new_s);
  // Attach it to the end of the linked list
  current_s=render->scenes;
  last_s   =current_s;
  while(current_s!=NULL){
    last_s   =current_s;
    current_s=current_s->next;
  }
  // Set first scene pointer
  if(last_s==NULL){
    render->first_scene=new_s;
    render->scenes     =new_s;
  }
  // Set last (ie last added) scene pointers
  //   (Carry last scene's last perspective as new starting perspective)
  else{
    copy_perspective(last_s->last_perspective,new_s->first_perspective);
    last_s->next=new_s;
  }
  render->last_scene=new_s;
  SID_log("Done.",SID_LOG_CLOSE);
}

