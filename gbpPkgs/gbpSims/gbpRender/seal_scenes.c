#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

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

