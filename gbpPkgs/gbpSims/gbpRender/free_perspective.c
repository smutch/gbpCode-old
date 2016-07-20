#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void free_perspective(perspective_info **perspective){
  perspective_info *current;
  perspective_info *next;
  SID_log("Freeing perspective...",SID_LOG_OPEN);
  // If this perspective has been allocated ...
  if((*perspective)!=NULL){
    // ... then free its linked list
    current=(*perspective);
    while(current!=NULL){
      next=current->next;
      SID_free(SID_FARG current);
      current=next;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);
}

