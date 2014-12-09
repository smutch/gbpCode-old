#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void free_render(render_info **render){
  int i_snap;
  SID_log("Freeing render structure...",SID_LOG_OPEN);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  free_camera(&((*render)->camera));
  free_scenes(&((*render)->scenes));
  if((*render)->plist_list!=NULL){
     for(i_snap=0;i_snap<(*render)->n_interpolate;i_snap++)
        free_plist((*render)->plist_list[i_snap]);
     SID_free(SID_FARG (*render)->plist_list);
     SID_free(SID_FARG (*render)->snap_list);
  }
  if((*render)->snap_a_list!=NULL)
     SID_free(SID_FARG (*render)->snap_a_list);
  SID_free(SID_FARG (*render)->kernel_radius);
  SID_free(SID_FARG (*render)->kernel_table);
  SID_free(SID_FARG (*render)->kernel_table_3d);
  free_mark_arguments(&((*render)->mark_arg_first));
  SID_free(SID_FARG (*render));
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);
}

