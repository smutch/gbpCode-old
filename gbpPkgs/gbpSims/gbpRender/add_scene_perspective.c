#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void add_scene_perspective(scene_info *scene){
  perspective_info *current_p;
  perspective_info *last_p;
  perspective_info *new_p;
  SID_log("Adding scene perspective...",SID_LOG_OPEN);
  // Create a new scene
  init_perspective(&new_p,RENDER_INIT_PERSPECTIVE);
  // Attach it to the end of the linked list
  current_p=scene->perspectives;
  last_p   =current_p;
  while(current_p!=NULL){
    last_p   =current_p;
    current_p=current_p->next;
  }
  // Set first perspective pointer
  if(last_p==NULL){
    scene->first_perspective=new_p;
    scene->perspectives     =new_p;
  }
  // Set last (ie last added) scene pointers; use previous perspective as new defaults
  else{
    copy_perspective(last_p,new_p);
    last_p->next=new_p;
  }
  scene->last_perspective=new_p;
  scene->n_perspectives++;
  SID_log("Done.",SID_LOG_CLOSE);
}

