#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpRender.h>

int main(int argc, char *argv[]){
  int          flag_continue;
  int          i_frame;
  int          i_frame_start;
  int          i_frame_stop;
  render_info *render=NULL;
  char         filename_script[MAX_FILENAME_LENGTH];
  int          mode;
  double       RGB_min;
  double       RGB_max;
  double       Y_min;
  double       Y_max;
  FILE        *fp_check;

  SID_init(&argc,&argv,NULL);

  // Parse cmd line arguments
  strcpy(filename_script,argv[1]);
  if(argc==4){
    i_frame_start=atoi(argv[2]);
    i_frame_stop =atoi(argv[3]);
  }
  else if(argc!=2)
    SID_trap_error("Invalid number of arguments.",ERROR_SYNTAX);

  SID_log("Rendering script file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_script);

  // Parse the script file and initialize the render structure
  parse_render_file(&render,filename_script);

  // Set image range (if not given on the cmd line)
  if(argc!=4){
    i_frame_start=0;
    i_frame_stop =render->n_frames-1;
  }

  // Decide if this is a fresh render based on whether the
  //   output directory exists or not
  if((fp_check=fopen((*render).filename_out_dir,"r"))!=NULL){
    mode=SET_RENDER_RESCALE;
    fclose(fp_check);
  }
  else
    mode=SET_RENDER_DEFAULT;
  
  // Loop over all the frames
  for(i_frame=i_frame_start;i_frame<=i_frame_stop;i_frame++){

    SID_log("Generating frame %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_frame+1,render->n_frames);

    // Set the render state
    set_render_state(render,i_frame,mode);
    
    // Render/read frame
    if(mode==SET_RENDER_RESCALE)
      read_frame(render,i_frame);
    else
      render_frame(render);

    // Write output
    if(SID.I_am_Master){
      if(mode==SET_RENDER_RESCALE)
        write_frame(render,i_frame,WRITE_IMAGE_DEFAULT|WRITE_IMAGE_PNG_ONLY);
      else
        write_frame(render,i_frame,WRITE_IMAGE_DEFAULT);
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
  
  // Clean-up 
  free_render(&render);
  
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
