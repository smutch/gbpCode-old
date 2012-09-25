#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void init_perspective(perspective_info **perspective,int mode){
  SID_log("Initializing perspective...",SID_LOG_OPEN);
  (*perspective)=(perspective_info *)SID_malloc(sizeof(perspective_info));
  (*perspective)->p_o[0]= 0.;
  (*perspective)->p_o[1]= 0.;
  (*perspective)->p_o[2]= 0.;
  (*perspective)->theta = 0.;
  (*perspective)->zeta  = 0.;
  (*perspective)->FOV   = 0.;
  if(check_mode_for_flag(mode,RENDER_INIT_EVOLVE)){
    (*perspective)->radius= 0.;
    (*perspective)->phi   = 0.;
    (*perspective)->time  = 0.;    
  }
  else{
    (*perspective)->radius= 1.;
    (*perspective)->phi   = 1.;
    (*perspective)->time  = 1.;
  }
  (*perspective)->next  = NULL;
  SID_log("Done.",SID_LOG_CLOSE);
}

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

void copy_perspective(perspective_info *from,perspective_info *to){
  SID_log("Copying perspective...",SID_LOG_OPEN);
  to->p_o[0]=from->p_o[0];  
  to->p_o[1]=from->p_o[1];  
  to->p_o[2]=from->p_o[2];  
  to->radius=from->radius;     
  to->FOV   =from->FOV;     
  to->theta =from->theta;   
  to->zeta  =from->zeta;    
  to->phi   =from->phi;     
  to->time  =from->time;    
  SID_log("Done.",SID_LOG_CLOSE);
}

void init_perspective_interp(perspective_interp_info **perspective_interp){
  SID_log("Initializing perspective interpolation structures...",SID_LOG_OPEN);
  (*perspective_interp)=(perspective_interp_info *)SID_malloc(sizeof(perspective_interp_info));
  (*perspective_interp)->p_o[0]=NULL;
  (*perspective_interp)->p_o[1]=NULL;
  (*perspective_interp)->p_o[2]=NULL;
  (*perspective_interp)->radius=NULL;
  (*perspective_interp)->FOV   =NULL;
  (*perspective_interp)->theta =NULL;
  (*perspective_interp)->zeta  =NULL;
  (*perspective_interp)->phi   =NULL;
  (*perspective_interp)->time  =NULL;
  SID_log("Done.",SID_LOG_CLOSE);
}

void free_perspective_interp(perspective_interp_info **perspective_interp){
  SID_log("Freeing perspective interpolation structures...",SID_LOG_OPEN);
  free_interpolate(SID_FARG (*perspective_interp)->p_o[0]);
  free_interpolate(SID_FARG (*perspective_interp)->p_o[1]);
  free_interpolate(SID_FARG (*perspective_interp)->p_o[2]);
  free_interpolate(SID_FARG (*perspective_interp)->theta);
  free_interpolate(SID_FARG (*perspective_interp)->FOV);
  free_interpolate(SID_FARG (*perspective_interp)->FOV);
  free_interpolate(SID_FARG (*perspective_interp)->radius);
  free_interpolate(SID_FARG (*perspective_interp)->phi);
  free_interpolate(SID_FARG (*perspective_interp)->time);
  SID_free(SID_FARG (*perspective_interp));
  SID_log("Done.",SID_LOG_CLOSE);
}

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

void seal_scenes(scene_info *scenes){
  scene_info       *current_scene;
  perspective_info *current_perspective;
  perspective_info *next;
  perspective_info *start;
  perspective_info *stop;
  perspective_info *last;
  double           *x;
  double           *y;
  int               i_frame;
  int               first_frame;
  int               last_frame;
  const gsl_interp_type  *interp_type;

  SID_log("Sealing scenes...",SID_LOG_OPEN);

  // Loop over the linked list of scenes ...
  current_scene=scenes;
  last_frame   =-1;
  while(current_scene!=NULL){
    first_frame=last_frame+1;
    last_frame =first_frame+current_scene->n_frames-1;
    
    if(!current_scene->sealed){
      current_scene->first_frame=first_frame;
      current_scene->last_frame =last_frame;
    
      // Create a final perspective if the only defined state is the initial state
      if(current_scene->n_perspectives<2)
        add_scene_perspective(current_scene);

      // Set interpolation type
      if(current_scene->n_perspectives<=2)
        interp_type=gsl_interp_linear;
      else
        interp_type=gsl_interp_cspline;

      // Add evolution to the perspectives (this way, last states get passed forward properly to the next scene)
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        current_perspective->p_o[0]+=current_scene->evolve->p_o[0]*(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->p_o[1]+=current_scene->evolve->p_o[1]*(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->p_o[2]+=current_scene->evolve->p_o[2]*(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->radius+=current_scene->evolve->radius*(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->FOV   +=current_scene->evolve->FOV   *(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->theta +=current_scene->evolve->theta *(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->zeta  +=current_scene->evolve->zeta  *(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->phi   +=current_scene->evolve->phi   *(double)(i_frame)/(double)(scenes->n_perspectives-1);
        current_perspective->time  +=current_scene->evolve->time  *(double)(i_frame)/(double)(scenes->n_perspectives-1);
        i_frame++;
        current_perspective=current_perspective->next;
      }

      // Create x-coordinate for frame interpolation
      x=(double *)SID_malloc(sizeof(double)*current_scene->n_perspectives);
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        x[i_frame]=(double)first_frame+(double)(i_frame*current_scene->n_frames)/(double)(current_scene->n_perspectives-1);
        current_perspective=current_perspective->next;
        i_frame++;
      }
      
      // Create y-coordinates for frame interpolation
      y=(double *)SID_malloc(sizeof(double)*current_scene->n_perspectives);
      // p_o
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->p_o[0];
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->p_o[0]));
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->p_o[1];
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->p_o[1]));
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->p_o[2];
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->p_o[2]));
      // theta
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->theta;
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->theta));
      // zeta
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->zeta;
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->zeta));
      // phi
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->phi;
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->phi));
      // Radius
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->radius;
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->radius));
      // FOV
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->FOV;
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->FOV));
      // time
      current_perspective=current_scene->perspectives;
      i_frame=0;
      while(current_perspective!=NULL){
        y[i_frame++]=current_perspective->time;
        current_perspective=current_perspective->next;
      }
      init_interpolate(x,y,current_scene->n_perspectives,interp_type,&(current_scene->interp->time));
      current_scene->sealed=TRUE;

      // Clean-up
      SID_free(SID_FARG x);
      SID_free(SID_FARG y);
    }
    current_scene=current_scene->next;
  } 
  SID_log("Done.",SID_LOG_CLOSE);
}

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

void init_camera(camera_info **camera, int mode){
  SID_log("Initializing camera...",SID_LOG_OPEN);

  // Initialize camera
  (*camera)=(camera_info *)SID_malloc(sizeof(camera_info));
  (*camera)->camera_mode=mode;
  
  // Initialize the perspective information for this camera
  init_perspective(&((*camera)->perspective),RENDER_INIT_PERSPECTIVE);

  // Initialze image information
  (*camera)->stereo_ratio =0.;
  (*camera)->RGB_mode     =0;
  (*camera)->RGB_gamma    =NULL;
  (*camera)->RGB_transfer =NULL;
  (*camera)->Y_mode       =0;
  (*camera)->Y_gamma      =NULL;
  (*camera)->Y_transfer   =NULL;
  (*camera)->Z_mode       =0;
  (*camera)->Z_gamma      =NULL;
  (*camera)->Z_transfer   =NULL;
  strcpy((*camera)->RGB_param,"");
  strcpy((*camera)->Y_param,  "");

  // Initialize image buffers
  (*camera)->colour_table    =4;
  (*camera)->image_RGB       =NULL;
  (*camera)->image_Y         =NULL;
  (*camera)->image_Z         =NULL;
  (*camera)->image_RGBY      =NULL;
  (*camera)->image_RGB_left  =NULL;
  (*camera)->image_Y_left    =NULL;
  (*camera)->image_Z_left    =NULL;
  (*camera)->image_RGBY_left =NULL;    
  (*camera)->image_RGB_right =NULL;
  (*camera)->image_Y_right   =NULL;
  (*camera)->image_Z_right   =NULL;
  (*camera)->image_RGBY_right=NULL;

  SID_log("Done.",SID_LOG_CLOSE);
}

void seal_render_camera(render_info *render){
  
  SID_log("Sealing camera...",SID_LOG_OPEN);

  // Decide if we're producing a stereo image based on stereo_factor
  if(render->camera->stereo_ratio>0.) 
    render->camera->camera_mode|=CAMERA_STEREO;

  // Initialize the perspective information for this camera
  copy_perspective(render->first_scene->first_perspective,render->camera->perspective);

  // Initialize image buffers (use colour_table=1 ... ie B&W ... for Y and Z images)
  init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGB));
  init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Y));
  if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL))
    init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Z));
  init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGBY));
  render->camera->mask_RGB =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
  render->camera->mask_Y   =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
  render->camera->mask_RGBY=(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
  if(check_mode_for_flag(render->camera->camera_mode,CAMERA_STEREO)){
    init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGB_left));
    init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Y_left));
    if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL))
      init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Z_left));
    init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGBY_left));    
    init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGB_right));
    init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Y_right));
    if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL))
      init_image(render->camera->width,render->camera->height,1,&(render->camera->image_Z_right));
    init_image(render->camera->width,render->camera->height,render->camera->colour_table,&(render->camera->image_RGBY_right));    
    render->camera->mask_RGB_left  =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_Y_left    =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_RGBY_left =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_RGB_right =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_Y_right   =(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
    render->camera->mask_RGBY_right=(char *)SID_malloc(sizeof(char)*render->camera->width*render->camera->height);
  }

  // Convert camera depth-range to Mpc/h
  render->camera->Z_range[0]*=M_PER_MPC/render->h_Hubble;
  render->camera->Z_range[1]*=M_PER_MPC/render->h_Hubble;
  SID_log("Done.",SID_LOG_CLOSE);
}

void free_camera(camera_info **camera){
  SID_log("Freeing camera...",SID_LOG_OPEN);
  free_perspective(&((*camera)->perspective));
  free_image(&((*camera)->image_RGB));
  free_image(&((*camera)->image_Y));
  free_image(&((*camera)->image_Z));
  free_image(&((*camera)->image_RGBY));
  SID_free(SID_FARG (*camera)->mask_RGB);
  SID_free(SID_FARG (*camera)->mask_Y);
  SID_free(SID_FARG (*camera)->mask_RGBY);
  if(check_mode_for_flag((*camera)->camera_mode,CAMERA_STEREO)){
    free_image(&((*camera)->image_RGB_left));
    free_image(&((*camera)->image_Y_left));
    free_image(&((*camera)->image_Z_left));
    free_image(&((*camera)->image_RGBY_left));
    free_image(&((*camera)->image_RGB_right));
    free_image(&((*camera)->image_Y_right));
    free_image(&((*camera)->image_Z_right));
    free_image(&((*camera)->image_RGBY_right));
    SID_free(SID_FARG (*camera)->mask_RGB_left);
    SID_free(SID_FARG (*camera)->mask_Y_left);
    SID_free(SID_FARG (*camera)->mask_RGBY_left);
    SID_free(SID_FARG (*camera)->mask_RGB_right);
    SID_free(SID_FARG (*camera)->mask_Y_right);
    SID_free(SID_FARG (*camera)->mask_RGBY_right);
  }
  if((*camera)->RGB_gamma!=NULL)
    free_interpolate(SID_FARG (*camera)->RGB_gamma);
  if((*camera)->RGB_transfer!=NULL)
    ADaPS_free(SID_FARG (*camera)->RGB_transfer);
  if((*camera)->Y_gamma!=NULL)
    free_interpolate(SID_FARG (*camera)->Y_gamma);
  if((*camera)->Y_transfer!=NULL)
    ADaPS_free(SID_FARG (*camera)->Y_transfer);
  if((*camera)->Z_gamma!=NULL)
    free_interpolate(SID_FARG (*camera)->Z_gamma);
  if((*camera)->Z_transfer!=NULL)
    ADaPS_free(SID_FARG (*camera)->Z_transfer);
  SID_free((void **)camera);
  SID_log("Done.",SID_LOG_CLOSE);
}

void init_render(render_info **render){
  int         flag_scatter;
  int         flag_no_velocities;
  SID_log("Initializing render structure...",SID_LOG_OPEN);

  // Allocate memory for rendering
  (*render)=(render_info *)SID_malloc(sizeof(render_info));

  // Initialize the camera
  init_camera(&((*render)->camera),CAMERA_DEFAULT);

  // Create the first scene
  init_scene(&((*render)->scenes));
  (*render)->first_scene=(*render)->scenes;
  (*render)->last_scene =(*render)->scenes;

  (*render)->n_frames           = 0;
  (*render)->snap_number_read   =-1;
  (*render)->snap_number        = 0;
  (*render)->n_snap_a_list      = 0;
  (*render)->snap_a_list        = NULL;
  (*render)->h_Hubble           = 1.;
  (*render)->near_field         = 0.;
  (*render)->kappa_absorption   =-1.;
  (*render)->kappa_transfer     = NULL;
  (*render)->flag_read_marked   = FALSE;
  (*render)->flag_comoving      = TRUE;
  (*render)->flag_force_periodic= FALSE;
  (*render)->flag_add_absorption= FALSE;
  (*render)->sealed             = FALSE;
  (*render)->mode               = MAKE_MAP_LOG;

  // Initialize the structure that holds the particle data
  init_plist(&((*render)->plist),NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  flag_scatter      =TRUE;
  flag_no_velocities=TRUE;
  ADaPS_store(&((*render)->plist.data),(void *)(&flag_scatter),      "flag_read_scatter", ADaPS_SCALAR_INT);
  ADaPS_store(&((*render)->plist.data),(void *)(&flag_no_velocities),"flag_no_velocities",ADaPS_SCALAR_INT);

  SID_log("Done.",SID_LOG_CLOSE);
}

void free_render(render_info **render){
  SID_log("Freeing render structure...",SID_LOG_OPEN);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  free_camera(&((*render)->camera));
  free_scenes(&((*render)->scenes));
  free_plist(&((*render)->plist));
  if((*render)->snap_a_list!=NULL)
    SID_free(SID_FARG (*render)->snap_a_list);
  SID_free(SID_FARG (*render));
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);
}

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

// Return true if the new perspective is sucessfully set
int set_render_state(render_info *render,int frame,int mode){
  scene_info       *current_scene;
  perspective_info *perspective;
  int               start_frame;
  int               stop_frame=0;
  int               r_val=FALSE;
  int               i_snap,snap_best;
  double            snap_diff,snap_diff_best;

  if(!render->sealed)
    seal_render(render);
    
  perspective  =render->camera->perspective;
  current_scene=render->scenes;
  while(current_scene!=NULL && r_val==FALSE){    
    r_val=TRUE;
    if(frame==current_scene->first_frame)
      copy_perspective(current_scene->first_perspective,perspective);
    else if(frame==current_scene->last_frame)
      copy_perspective(current_scene->last_perspective, perspective);
    else if(frame>=current_scene->first_frame && frame<current_scene->last_frame){
      perspective->p_o[0]=interpolate(current_scene->interp->p_o[0],(double)frame);
      perspective->p_o[1]=interpolate(current_scene->interp->p_o[1],(double)frame);
      perspective->p_o[2]=interpolate(current_scene->interp->p_o[2],(double)frame);
      perspective->time  =interpolate(current_scene->interp->time,  (double)frame);
      perspective->radius=interpolate(current_scene->interp->radius,(double)frame);
      perspective->FOV   =interpolate(current_scene->interp->FOV,   (double)frame);
      perspective->theta =interpolate(current_scene->interp->theta, (double)frame);
      perspective->zeta  =interpolate(current_scene->interp->zeta,  (double)frame);
      perspective->phi   =interpolate(current_scene->interp->phi,   (double)frame);
    }
    else
      r_val=FALSE;
    current_scene=current_scene->next;
  }
  if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL)) 
    perspective->radius*=perspective->phi;
  else
    perspective->radius=1e8*M_PER_MPC;

  perspective->FOV   *=perspective->phi;
  perspective->p_c[0] =perspective->p_o[0]+perspective->radius*cos(perspective->zeta)*sin(perspective->theta);
  perspective->p_c[1] =perspective->p_o[1]+perspective->radius*cos(perspective->zeta)*cos(perspective->theta);
  perspective->p_c[2] =perspective->p_o[2]+perspective->radius*sin(perspective->zeta);
  perspective->d_o    =sqrt(pow(perspective->p_o[0]-perspective->p_c[0],2)+
                            pow(perspective->p_o[1]-perspective->p_c[1],2)+
                            pow(perspective->p_o[2]-perspective->p_c[2],2));
 
  // Perform snapshot and smooth-file reading
  if(!check_mode_for_flag(mode,SET_RENDER_RESCALE)){
    // Determine which snapshot(s) to use
    if(render->snap_a_list!=NULL && render->n_snap_a_list>0){
      SID_log("Selecting snapshot for t=%lf...",SID_LOG_OPEN,perspective->time);
      snap_best     =0;
      snap_diff_best=1e60;
      for(i_snap=0;i_snap<render->n_snap_a_list;i_snap++){
        snap_diff=fabs(perspective->time-render->snap_a_list[i_snap]);
        if(snap_diff<snap_diff_best){
          snap_diff_best=snap_diff;
          snap_best     =i_snap;
        }
      }
      render->snap_number=snap_best;
      SID_log("snap=%d is best with t=%lf...Done.",SID_LOG_CLOSE,render->snap_number,render->snap_a_list[render->snap_number]);
    }

    // Check (and read if necessary) snapshots and smooth files here
    if(render->snap_number_read!=render->snap_number){
      if(render->flag_read_marked){
         read_mark_file(&(render->plist),"mark",render->mark_filename_root,MARK_LIST_ONLY);
         render->flag_read_marked=FALSE;
      }
      if(render->flag_add_absorption){
         ADaPS_remove(&(render->plist.data),"flag_read_scatter");
         read_gadget_binary(render->snap_filename_root,render->snap_number,&(render->plist),READ_GADGET_DEFAULT);
      }
      else
         read_gadget_binary_render(render->snap_filename_root,render->snap_number,&(render->plist));
      render->h_Hubble=((double *)ADaPS_fetch(render->plist.data,"h_Hubble"))[0];
      read_smooth(&(render->plist),render->smooth_filename_root,render->snap_number,SMOOTH_DEFAULT);
      render->snap_number_read=render->snap_number;
    }
  }
  
  // Convert [Mpc/h] -> SI
  perspective->p_o[0]*=M_PER_MPC/render->h_Hubble;
  perspective->p_o[1]*=M_PER_MPC/render->h_Hubble;
  perspective->p_o[2]*=M_PER_MPC/render->h_Hubble;
  perspective->p_c[0]*=M_PER_MPC/render->h_Hubble;
  perspective->p_c[1]*=M_PER_MPC/render->h_Hubble;
  perspective->p_c[2]*=M_PER_MPC/render->h_Hubble;
  perspective->d_o   *=M_PER_MPC/render->h_Hubble;
  perspective->radius*=M_PER_MPC/render->h_Hubble;
  perspective->FOV   *=M_PER_MPC/render->h_Hubble;
  
  return(r_val);
}

void parse_render_file(render_info **render, char *filename){
  FILE   *fp;
  int     n_lines;
  int     i_line;
  char   *line=NULL;
  size_t  line_length=0;
  int     i_word,j_word;
  char    temp_word[1024];
  char    line_print[256];
  char    command[64];
  char    parameter[64];
  char    variable[64];
  char    image[64];
  char    c_value[64];
  double  d_value;
  int     i_value;
  int     i;
  int     flag;
  int     n_transfer;
  double      *transfer_array_x;
  double      *transfer_array_y;
  interp_info *temp_interp;
  FILE        *fp_list;
  
  SID_log("Parsing render file {%s}...",SID_LOG_OPEN,filename);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
  if((fp=fopen(filename,"r"))==NULL)
    SID_trap_error("Could not open render file {%s}",ERROR_IO_OPEN,filename);
  init_render(render);
  sprintf((*render)->filename_out_dir,"%s.images/",filename);
  n_lines=count_lines_data(fp);
  for(i_line=0;i_line<n_lines;i_line++){
    SID_log("Processing line %d of %d...",SID_LOG_OPEN,i_line+1,n_lines);SID_Barrier(SID.COMM_WORLD);
    grab_next_line_data(fp,&line,&line_length);
    i_word=1;
    grab_word(line,i_word++,command);
    if(strlen(line)>0){
      sprintf(line_print,"%s",line);
      // Interpret command
      if(!strcmp(command,"set")){
        grab_word(line,i_word++,parameter);
        if(!strcmp(parameter,"scene")){
          grab_word(line,i_word++,variable);
          if(!strcmp(variable,"n_frames")){
            grab_int(line,i_word++,&i_value);
            (*render)->last_scene->n_frames=i_value;
          }
          else if(!strcmp(variable,"FOV")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->FOV=d_value;
          }
          else if(!strcmp(variable,"p_o")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->p_o[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->p_o[1]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->p_o[2]=d_value;
          }
          else if(!strcmp(variable,"radius")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->radius=d_value;
          }
          else if(!strcmp(variable,"theta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->theta=d_value;
          }
          else if(!strcmp(variable,"phi")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->phi=d_value;
          }
          else if(!strcmp(variable,"zeta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->zeta=d_value;
          }
          else if(!strcmp(variable,"time")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->last_perspective->time=d_value;
          }
          else
            SID_trap_error("Unknown variable {%s} for parameter {%s} for command {%s} on line %d",ERROR_LOGIC,variable,parameter,command,i_line);
        }
        else if(!strcmp(parameter,"snap_file"))
          grab_word(line,i_word++,(*render)->snap_filename_root);
        else if(!strcmp(parameter,"mark_file")){
          grab_word(line,i_word++,(*render)->mark_filename_root);
          (*render)->flag_read_marked=TRUE;
        }
        else if(!strcmp(parameter,"smooth_file"))
          grab_word(line,i_word++,(*render)->smooth_filename_root);
        else if(!strcmp(parameter,"snap_a_list_file")){
          SID_log("Reading snapshot a_list...",SID_LOG_OPEN);
          grab_word(line,i_word++,(*render)->snap_a_list_filename);
          if((fp_list=fopen((*render)->snap_a_list_filename,"r"))==NULL)
            SID_trap_error("Could not open snapshot a_list {%s}",ERROR_IO_OPEN,(*render)->snap_a_list_filename);
          (*render)->n_snap_a_list=count_lines_data(fp_list);
          (*render)->snap_a_list  =(double *)SID_malloc(sizeof(double)*(*render)->n_snap_a_list);
          for(i=0;i<(*render)->n_snap_a_list;i++){
            grab_next_line_data(fp_list,&line,&line_length);
            grab_double(line,1,&((*render)->snap_a_list[i]));
          }
          fclose(fp_list);
          for(i=0;i<(*render)->n_snap_a_list;i++)
            SID_log("a[%3d]=%lf",SID_LOG_COMMENT,i,(*render)->snap_a_list[i]);
          SID_log("Done.",SID_LOG_CLOSE);
        }
        else if(!strcmp(parameter,"h_Hubble"))
          grab_double(line,i_word++,&((*render)->h_Hubble));
        else if(!strcmp(parameter,"near_field"))
          grab_double(line,i_word++,&((*render)->near_field));
        else if(!strcmp(parameter,"kappa_absorption")){
          if((*render)->flag_add_absorption)
             SID_trap_error("There are conflicting absorption criteria.",ERROR_LOGIC);
          (*render)->flag_add_absorption=TRUE;
          grab_double(line,i_word++,&((*render)->kappa_absorption));
        }
        else if(!strcmp(parameter,"kappa_transfer")){
          if((*render)->flag_add_absorption)
             SID_trap_error("There are conflicting absorption criteria.",ERROR_LOGIC);
          (*render)->flag_add_absorption=TRUE;
          i_word=set_transfer_function(line,i_word,&((*render)->kappa_transfer));
        }
        else if(!strcmp(parameter,"snap_number"))
          grab_int(line,i_word++,&((*render)->snap_number));
        else if(!strcmp(parameter,"flag_comoving"))
          grab_int(line,i_word++,&((*render)->flag_comoving));
        else if(!strcmp(parameter,"force_periodic"))
          (*render)->flag_force_periodic=TRUE;
        else if(!strcmp(parameter,"camera")){
          grab_word(line,i_word++,variable);
          if(!strcmp(variable,"size")){
            grab_int(line,i_word++,&i_value);
            (*render)->camera->width =i_value;
            grab_int(line,i_word++,&i_value);
            (*render)->camera->height=i_value;
          }
          else if(!strcmp(variable,"stereo_factor")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->stereo_ratio=d_value;
          }
          else if(!strcmp(variable,"colour_table")){
            grab_int(line,i_word++,&i_value);
            (*render)->camera->colour_table=i_value;
          }
          else if(!strcmp(variable,"plane_parallel")){
            (*render)->camera->camera_mode|=CAMERA_PLANE_PARALLEL|CAMERA_MONO;
          }
          else if(!strcmp(variable,"RGB_param")){
            grab_word(line,i_word++,c_value);
            strcpy((*render)->camera->RGB_param,c_value);
          }
          else if(!strcmp(variable,"RGB_range")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->RGB_range[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->camera->RGB_range[1]=d_value;
          }
          else if(!strcmp(variable,"RGB_gamma"))
            i_word=set_transfer_function(line,i_word,&((*render)->camera->RGB_gamma));
          else if(!strcmp(variable,"RGB_transfer")){
            grab_word(line,i_word++,parameter);
            grab_word(line,i_word++,temp_word);
            if(!strcmp(temp_word,"log") || !strcmp(temp_word,"LOG"))
              flag=TRUE;
            else if(!strcmp(temp_word,"linear") || !strcmp(temp_word,"LINEAR"))
              flag=FALSE;
            else
              SID_trap_error("log/linear flag not set to 'log' or 'linear' {%s}",ERROR_LOGIC,temp_word);
            i_word=set_transfer_function(line,i_word,&temp_interp);
            // Add to the transfer function list
            ADaPS_store(&((*render)->camera->RGB_transfer),(void *)temp_interp,parameter,ADaPS_DEFAULT);
            // If this is a log-defined transfer function, store that fact
            if(flag)
              ADaPS_store(&((*render)->camera->RGB_transfer),(void *)&flag,"%s_log",ADaPS_SCALAR_INT,parameter);
          }
          else if(!strcmp(variable,"Y_param")){
            grab_word(line,i_word++,c_value);
            strcpy((*render)->camera->Y_param,c_value);
          }
          else if(!strcmp(variable,"Y_range")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Y_range[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Y_range[1]=d_value;
          }
          else if(!strcmp(variable,"Y_gamma")){
            n_transfer=count_words(line)-i_word+1;
            if(n_transfer>2){
              transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
              transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
              if((*render)->camera->Y_gamma!=NULL)
                free_interpolate(SID_FARG (*render)->camera->Y_gamma);
              for(j_word=0;j_word<n_transfer;j_word++){
                grab_double(line,i_word++,&d_value);
                transfer_array_x[j_word]=255.*(double)j_word/(double)(n_transfer-1);
                transfer_array_y[j_word]=d_value;
              }
              init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,&((*render)->camera->Y_gamma));
              SID_free(SID_FARG transfer_array_x);
              SID_free(SID_FARG transfer_array_y);
            }
            else
              SID_log_warning("Gamma arrays bust be >2 elements long.",ERROR_LOGIC);
          }
          else if(!strcmp(variable,"Y_transfer")){
            grab_word(line,i_word++,parameter);
            grab_word(line,i_word++,temp_word);
            if(!strcmp(temp_word,"log") || !strcmp(temp_word,"LOG"))
              flag=TRUE;
            else if(!strcmp(temp_word,"linear") || !strcmp(temp_word,"LINEAR"))
              flag=FALSE;
            else
              SID_trap_error("log/linear flag not set to 'log' or 'linear' {%s}",ERROR_LOGIC,temp_word);
            n_transfer      =count_words(line)-i_word+1;
            SID_log("parameter ={%s}",SID_LOG_COMMENT,parameter);
            SID_log("log/lin   ={%s}",SID_LOG_COMMENT,temp_word);
            SID_log("n_transfer=%d",SID_LOG_COMMENT,n_transfer);
            if(n_transfer>2){
              n_transfer+=2; // We need to add low/hi interpolation anchors
              transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
              transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
              for(j_word=1;j_word<n_transfer-1;j_word++){
                grab_word(line,i_word++,temp_word);
                search_and_replace(temp_word,","," ");
                if(count_words(temp_word)!=2)
                  SID_trap_error("Error in formatting of transfer array {%s}{%s}",ERROR_LOGIC,line,temp_word);
                grab_double(temp_word,1,&(transfer_array_x[j_word]));
                grab_double(temp_word,2,&(transfer_array_y[j_word]));
                SID_log("%le %le",SID_LOG_COMMENT,transfer_array_x[j_word],transfer_array_y[j_word]);
              }
              // Create low/hi interpolation anchors
              transfer_array_y[0]           =transfer_array_y[1];
              transfer_array_y[n_transfer-1]=transfer_array_y[n_transfer-2];
              if(transfer_array_x[1]<0.)
                transfer_array_x[0]=10.*transfer_array_x[1];
              else
                transfer_array_x[0]=0.1*transfer_array_x[1];
              if(transfer_array_x[n_transfer-2]<0.)
                transfer_array_x[n_transfer-1]= 0.1*transfer_array_x[n_transfer-2];
              else
                transfer_array_x[n_transfer-1]=10.0*transfer_array_x[n_transfer-2];
              // Create interpolation array
              init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,&temp_interp);
              // Add to the transfer function list
              ADaPS_store(&((*render)->camera->Y_transfer),(void *)temp_interp,parameter,ADaPS_DEFAULT);
              // If this is a log-defined transfer function, store that fact
              if(flag)
                ADaPS_store(&((*render)->camera->Y_transfer),(void *)&flag,"%s_log",ADaPS_SCALAR_INT,parameter);
              // Free temporary arrays
              SID_free(SID_FARG transfer_array_x);
              SID_free(SID_FARG transfer_array_y);
            }
            else
              SID_log_warning("Transfer arrays must be >2 elements long.",ERROR_LOGIC);
          }
          else if(!strcmp(variable,"Z_range")){
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Z_range[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->camera->Z_range[1]=d_value;
          }
          else if(!strcmp(variable,"Z_gamma")){
            n_transfer=count_words(line)-i_word+1;
            if(n_transfer>2){
              transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
              transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
              if((*render)->camera->Z_gamma!=NULL)
                free_interpolate(SID_FARG (*render)->camera->Z_gamma);
              for(j_word=0;j_word<n_transfer;j_word++){
                grab_double(line,i_word++,&d_value);
                transfer_array_x[j_word]=255.*(double)j_word/(double)(n_transfer-1);
                transfer_array_y[j_word]=d_value;
              }
              init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,&((*render)->camera->Z_gamma));
              SID_free(SID_FARG transfer_array_x);
              SID_free(SID_FARG transfer_array_y);
            }
            else
              SID_log_warning("Gamma arrays bust be >2 elements long.",ERROR_LOGIC);
          }
          else if(!strcmp(variable,"Z_transfer")){
            grab_word(line,i_word++,parameter);
            grab_word(line,i_word++,temp_word);
            if(!strcmp(temp_word,"log") || !strcmp(temp_word,"LOG"))
              flag=TRUE;
            else if(!strcmp(temp_word,"linear") || !strcmp(temp_word,"LINEAR"))
              flag=FALSE;
            else
              SID_trap_error("log/linear flag not set to 'log' or 'linear' {%s}",ERROR_LOGIC,temp_word);
            n_transfer      =count_words(line)-i_word+1;
            SID_log("parameter ={%s}",SID_LOG_COMMENT,parameter);
            SID_log("log/lin   ={%s}",SID_LOG_COMMENT,temp_word);
            SID_log("n_transfer=%d",SID_LOG_COMMENT,n_transfer);
            if(n_transfer>2){
              n_transfer+=2; // We need to add low/hi interpolation anchors
              transfer_array_x=(double *)SID_malloc(sizeof(double)*n_transfer);
              transfer_array_y=(double *)SID_malloc(sizeof(double)*n_transfer);
              for(j_word=1;j_word<n_transfer-1;j_word++){
                grab_word(line,i_word++,temp_word);
                search_and_replace(temp_word,","," ");
                if(count_words(temp_word)!=2)
                  SID_trap_error("Error in formatting of transfer array {%s}{%s}",ERROR_LOGIC,line,temp_word);
                grab_double(temp_word,1,&(transfer_array_x[j_word]));
                grab_double(temp_word,2,&(transfer_array_y[j_word]));
                SID_log("%le %le",SID_LOG_COMMENT,transfer_array_x[j_word],transfer_array_y[j_word]);
              }
              // Create low/hi interpolation anchors
              transfer_array_y[0]           =transfer_array_y[1];
              transfer_array_y[n_transfer-1]=transfer_array_y[n_transfer-2];
              if(transfer_array_x[1]<0.)
                transfer_array_x[0]=10.*transfer_array_x[1];
              else
                transfer_array_x[0]=0.1*transfer_array_x[1];
              if(transfer_array_x[n_transfer-2]<0.)
                transfer_array_x[n_transfer-1]= 0.1*transfer_array_x[n_transfer-2];
              else
                transfer_array_x[n_transfer-1]=10.0*transfer_array_x[n_transfer-2];
              // Create interpolation array
              init_interpolate(transfer_array_x,transfer_array_y,n_transfer,gsl_interp_cspline,&temp_interp);
              // Add to the transfer function list
              ADaPS_store(&((*render)->camera->Z_transfer),(void *)temp_interp,parameter,ADaPS_DEFAULT);
              // If this is a log-defined transfer function, store that fact
              if(flag)
                ADaPS_store(&((*render)->camera->Z_transfer),(void *)&flag,"%s_log",ADaPS_SCALAR_INT,parameter);
              // Free temporary arrays
              SID_free(SID_FARG transfer_array_x);
              SID_free(SID_FARG transfer_array_y);
            }
            else
              SID_log_warning("Transfer arrays must be >2 elements long.",ERROR_LOGIC);
          }
          else
            SID_trap_error("Unknown variable {%s} for parameter {%s} for command {%s} on line %d",ERROR_LOGIC,variable,parameter,command,i_line);
        }
        else if(!strcmp(parameter,"evolve")){
          grab_word(line,i_word++,variable);
          if(!strcmp(variable,"FOV")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->FOV=d_value;
          }
          else if(!strcmp(variable,"p_o")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->p_o[0]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->p_o[1]=d_value;
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->p_o[2]=d_value;
          }
          else if(!strcmp(variable,"radius")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->radius=d_value;
          }
          else if(!strcmp(variable,"theta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->theta=d_value;
          }
          else if(!strcmp(variable,"phi")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->phi=d_value;
          }
          else if(!strcmp(variable,"zeta")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->zeta=d_value;
          }
          else if(!strcmp(variable,"time")){
            grab_double(line,i_word++,&d_value);
            (*render)->last_scene->evolve->time=d_value;
          }
          else
            SID_trap_error("Unknown variable {%s} for parameter {%s} for command {%s} on line %d",ERROR_LOGIC,variable,parameter,command,i_line);
        }      
        else
          SID_trap_error("Unknown parameter {%s} for command {%s} on line %d",ERROR_LOGIC,parameter,command,i_line);
      }
      else if(!strcmp(command,"add")){
        grab_word(line,i_word++,parameter);
        if(!strcmp(parameter,"scene"))
          add_render_scene((*render));
        else if(!strcmp(parameter,"perspective"))
          add_scene_perspective((*render)->last_scene);
        else
          SID_trap_error("Unknown parameter {%s} for command {%s} on line %d",ERROR_LOGIC,parameter,command,i_line);
      }
      else
        SID_trap_error("Unknown render command {%s} on line %d",ERROR_LOGIC,command,i_line);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  fclose(fp);
  if(!(*render)->sealed)
    seal_render(*render);
  if(line!=NULL)
     SID_free(SID_FARG line);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);
}

void read_frame(render_info *render,int frame){
  char filename_RGB[256];
  char filename_Y[256];
  char filename_Z[256];
  char filename_RGBY[256];

  SID_log("Writing rendered frame...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Create directory if needed
  mkdir(render->filename_out_dir,02755);
  
  sprintf(filename_RGB, "%s/RGB_M_%05d", render->filename_out_dir,frame);
  sprintf(filename_Y,   "%s/Y_M_%05d",   render->filename_out_dir,frame);
  sprintf(filename_Z,   "%s/Z_M_%05d",   render->filename_out_dir,frame);
  sprintf(filename_RGBY,"%s/RGBY_M_%05d",render->filename_out_dir,frame);
  read_image(render->camera->image_RGB, filename_RGB);
  read_image(render->camera->image_Y,   filename_Y);
  read_image(render->camera->image_Z,   filename_Z);
  read_image(render->camera->image_RGBY,filename_RGBY);

  // Read stereo-images
  if(check_mode_for_flag(render->camera->camera_mode,CAMERA_STEREO)){
    // Left
    sprintf(filename_RGB, "%s/RGB_L_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_L_%05d",render->filename_out_dir,frame);
    read_image(render->camera->image_RGB_left,  filename_RGB);
    read_image(render->camera->image_Y_left,    filename_Y);
    read_image(render->camera->image_Z_left,    filename_Z);
    read_image(render->camera->image_RGBY_left, filename_RGBY);

    // Right
    sprintf(filename_RGB, "%s/RGB_R_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_R_%05d",render->filename_out_dir,frame);
    read_image(render->camera->image_RGB_right, filename_RGB);
    read_image(render->camera->image_Y_right,   filename_Y);
    read_image(render->camera->image_Z_right,   filename_Z);
    read_image(render->camera->image_RGBY_right,filename_RGBY);
  }  

  set_frame(render->camera);

  SID_log("Done.",SID_LOG_CLOSE);
}

void write_frame(render_info *render,int frame,int mode){
  char filename_RGB[256];
  char filename_Y[256];
  char filename_Z[256];
  char filename_RGBY[256];

  SID_log("Writing rendered frame...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Create directory if needed
  mkdir(render->filename_out_dir,02755);

  // Write mono-images
  set_frame(render->camera);
  
  sprintf(filename_RGB, "%s/RGB_M_%05d", render->filename_out_dir,frame);
  sprintf(filename_Y,   "%s/Y_M_%05d",   render->filename_out_dir,frame);
  sprintf(filename_Z,   "%s/Z_M_%05d",   render->filename_out_dir,frame);
  sprintf(filename_RGBY,"%s/RGBY_M_%05d",render->filename_out_dir,frame);
  write_image(render->camera->image_RGB, filename_RGB, mode);
  write_image(render->camera->image_Y,   filename_Y,   mode);
  if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL))
     write_image(render->camera->image_Z,   filename_Z,   mode);
  write_image(render->camera->image_RGBY,filename_RGBY,mode);

  // Write stereo-images
  if(check_mode_for_flag(render->camera->camera_mode,CAMERA_STEREO)){
    // Left
    sprintf(filename_RGB, "%s/RGB_L_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_L_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_L_%05d",render->filename_out_dir,frame);
    write_image(render->camera->image_RGB_left,  filename_RGB, mode);
    write_image(render->camera->image_Y_left,    filename_Y,   mode);
    if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL))
      write_image(render->camera->image_Z_left,    filename_Z,   mode);
    write_image(render->camera->image_RGBY_left, filename_RGBY,mode);

    // Right
    sprintf(filename_RGB, "%s/RGB_R_%05d", render->filename_out_dir,frame);
    sprintf(filename_Y,   "%s/Y_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_Z,   "%s/Z_R_%05d",   render->filename_out_dir,frame);
    sprintf(filename_RGBY,"%s/RGBY_R_%05d",render->filename_out_dir,frame);
    write_image(render->camera->image_RGB_right, filename_RGB, mode);
    write_image(render->camera->image_Y_right,   filename_Y,   mode);
    if(!check_mode_for_flag(render->camera->camera_mode,CAMERA_PLANE_PARALLEL))
      write_image(render->camera->image_Z_right,   filename_Z,   mode);
    write_image(render->camera->image_RGBY_right,filename_RGBY,mode);
  }  
  SID_log("Done.",SID_LOG_CLOSE);
}

void set_frame(camera_info *camera){
  int         i_x,i_y,i_pixel;
  int         i_image;
  double      image_min,image_max;
  double      image_range;
  char        pixel_value;
  double      brightness;
  image_info *image_RGB;
  image_info *image_Y;
  image_info *image_Z;
  image_info *image_RGBY;
  image_info *image;
  double     *values;

  // Loop over each set of images
  for(i_image=0;i_image<3;i_image++){
    switch(i_image){
    case 0:
      image_RGB =camera->image_RGB;      
      image_Y   =camera->image_Y;      
      image_Z   =camera->image_Z;      
      image_RGBY=camera->image_RGBY;
      break;
    case 1:      
      image_RGB =camera->image_RGB_left;      
      image_Y   =camera->image_Y_left;      
      image_Z   =camera->image_Z_left;      
      image_RGBY=camera->image_RGBY_left;
      break;
    case 2:
      image_RGB =camera->image_RGB_right;      
      image_Y   =camera->image_Y_right;
      image_Z   =camera->image_Z_right;
      image_RGBY=camera->image_RGBY_right;
      break;
    }

    // Set RGB Image
    if(image_RGB!=NULL)
       set_image_RGB(image_RGB,camera->RGB_range[0],camera->RGB_range[1]);
  
    // Set Y Image
    if(image_Y!=NULL)
       set_image_RGB(image_Y,camera->Y_range[0],camera->Y_range[1]);

    // Set Z Image
    if(image_Z!=NULL)
       set_image_RGB(image_Z,camera->Z_range[0],camera->Z_range[1]);

    // Set RGBY Image
    if(image_RGBY!=NULL && image_RGB!=NULL && image_Y!=NULL)
       set_image_RGBY(image_RGBY,
                      image_RGB,
                      image_Y,
                      camera->RGB_range[0],
                      camera->RGB_range[1],
                      camera->Y_range[0],
                      camera->Y_range[1]);

  }
}

void set_render_scale(render_info *render,double RGB_min,double RGB_max,double Y_min,double Y_max,double Z_min,double Z_max){
  render->camera->RGB_range[0]=RGB_min;
  render->camera->RGB_range[1]=RGB_max;
  render->camera->Y_range[0]  =Y_min;
  render->camera->Y_range[1]  =Y_max;
  render->camera->Z_range[0]  =Z_min;
  render->camera->Z_range[1]  =Z_max;
}
