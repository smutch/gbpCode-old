#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void free_scenes(scene_info **scene){
  scene_info *current;
  scene_info *next;
  SID_log("Freeing scene...",SID_LOG_OPEN);
  // If this scene has been allocated ...
  if((*scene)!=NULL){
    // ... free its linked list
    current=(*scene);
    while(current!=NULL){
      next=current->next;
      free_perspective(&(current->perspectives));
      free_perspective(&(current->evolve));
      free_perspective_interp(&(current->interp));
      SID_free(SID_FARG current);
      current=next;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);
}

