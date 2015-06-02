#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

// Return true if the new perspective is sucessfully set
int set_render_state(render_info *render,int frame,int mode){
  scene_info       *current_scene;
  perspective_info *perspective;
  int               start_frame;
  int               stop_frame=0;
  int               r_val=FALSE;
  int               i_snap,j_snap,snap_best;
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
    perspective->radius=1e8;

  perspective->FOV   *=perspective->phi;
  perspective->p_c[0] =perspective->p_o[0]+perspective->radius*cos(perspective->zeta)*sin(perspective->theta);
  perspective->p_c[1] =perspective->p_o[1]+perspective->radius*cos(perspective->zeta)*cos(perspective->theta);
  perspective->p_c[2] =perspective->p_o[2]+perspective->radius*sin(perspective->zeta);
  perspective->d_o    =sqrt(pow(perspective->p_o[0]-perspective->p_c[0],2)+
                            pow(perspective->p_o[1]-perspective->p_c[1],2)+
                            pow(perspective->p_o[2]-perspective->p_c[2],2));

  // Perform snapshot and smooth-file reading
  if(!check_mode_for_flag(mode,SET_RENDER_RESCALE)){
    int *snap_list;
    snap_list=(int *)SID_malloc(sizeof(int)*render->n_interpolate);
    // Determine which snapshot(s) to use
    if(render->snap_a_list!=NULL && render->n_snap_a_list>0){
      if(render->n_interpolate>1)
         SID_log("Selecting snapshots for t=%lf...",SID_LOG_OPEN,perspective->time);
      else
         SID_log("Selecting snapshot for t=%lf...",SID_LOG_OPEN,perspective->time);
      pick_best_snap(perspective->time,render->snap_a_list,render->n_snap_a_list,&snap_best,&snap_diff_best);
      if(render->n_interpolate>1){
         // Determine the list of best snapshots to use
         if(snap_diff_best<0)
            snap_list[0]=MAX(0,snap_best-render->n_interpolate/2);
         else
            snap_list[0]=MAX(0,1+snap_best-render->n_interpolate/2);
         for(i_snap=1;i_snap<render->n_interpolate;i_snap++)
            snap_list[i_snap]=snap_list[i_snap-1]+1;
         if(snap_list[render->n_interpolate-1]>=render->n_snap_a_list){
            int j_snap;
            for(i_snap=render->n_interpolate-1,j_snap=0;i_snap>=0;i_snap--,j_snap++)
               snap_list[i_snap]=render->n_snap_a_list-1-j_snap;
         }
         for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
            SID_log("%03d (time=%8.6f)",SID_LOG_COMMENT,snap_list[i_snap],render->snap_a_list[snap_list[i_snap]]);
         // Set interpolation factor (a_i-a_0)/(a_1-a_0)
         // WARNING: THIS ASSUMES THAT n_interpolate==2
         render->f_interpolate=(perspective->time-render->snap_a_list[snap_list[0]])/
                               (render->snap_a_list[snap_list[1]]-render->snap_a_list[snap_list[0]]);
         SID_log("f_interpolate=%le",SID_LOG_COMMENT,render->f_interpolate);
         if(render->n_interpolate!=2)
            SID_trap_error("n_interpolate>2 not supported (yet).",ERROR_NONE);

         SID_log("Done.",SID_LOG_CLOSE);
      }
      else if(render->n_interpolate<=0)
         SID_trap_error("An invalid value for n_interpolate (%d) has been set.",ERROR_LOGIC,render->n_interpolate);
      else{
         snap_list[0]=snap_best;
         SID_log("snap=%d is best with t=%lf...Done.",SID_LOG_CLOSE,snap_list[0],render->snap_a_list[snap_best]);
      }
    }

    // Initialize the array of pointers used to hold the snapshots
    if(render->plist_list==NULL){
       render->plist_list=(plist_info **)SID_malloc(sizeof(plist_info *)*(render->n_interpolate));
       render->snap_list =(int *)SID_calloc(sizeof(int)*(render->n_interpolate));
       for(i_snap=0;i_snap<render->n_interpolate;i_snap++){
          render->snap_list[i_snap] =-1;
          render->plist_list[i_snap]=(plist_info *)SID_malloc(sizeof(plist_info));
          init_plist(render->plist_list[i_snap],NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
       }
    }

    // Move currently loaded snapshots to new locations in the list if necessary
    for(i_snap=0;i_snap<render->n_interpolate;i_snap++){
       if(snap_list[i_snap]!=render->snap_list[i_snap]){
          // Free plists we are done with before loading a new one to avoid
          //   unnecessarily doubling-up on RAM usage
          if(render->plist_list[i_snap]!=NULL)
             ADaPS_free(SID_FARG render->plist_list[i_snap]->data);
          render->snap_list[i_snap]=-1;
          for(j_snap=i_snap+1;j_snap<render->n_interpolate;j_snap++){
             if(snap_list[i_snap]==render->snap_list[j_snap]){
                render->plist_list[i_snap]->data=render->plist_list[j_snap]->data;
                render->snap_list[i_snap]       =render->snap_list[j_snap];
                render->plist_list[j_snap]->data=NULL;
                render->snap_list[j_snap]       =-1;
             }
          }
       }
    }

    // Check and (if necessary) read snapshots and smooth files here
    for(i_snap=0;i_snap<render->n_interpolate;i_snap++){
       if(render->snap_list[i_snap]<0){
          render->snap_list[i_snap]=snap_list[i_snap];
          read_gadget_binary_render(render->snap_filename_root,render->snap_list[i_snap],render->plist_list[i_snap],READ_GADGET_RENDER_SCATTER);
          if(render->n_interpolate>1)
             read_smooth(render->plist_list[i_snap],render->smooth_filename_root,render->snap_list[i_snap],
                         SMOOTH_DEFAULT|READ_SMOOTH_LOG_SIGMA|READ_SMOOTH_LOG_RHO); // this is to speed-up logarythmic interpolation
          else
             read_smooth(render->plist_list[i_snap],render->smooth_filename_root,render->snap_list[i_snap],SMOOTH_DEFAULT);
       }
    }
    SID_free(SID_FARG snap_list);

    // Convert [Mpc/h] -> SI
    if(render->plist_list!=NULL){
       if(ADaPS_exist(render->plist_list[0]->data,"h_Hubble"))
          render->h_Hubble=((double *)ADaPS_fetch(render->plist_list[0]->data,"h_Hubble"))[0];
       else
          render->h_Hubble=1.;
    }
    else
       render->h_Hubble=1.;

    // Mark particles
    perform_marking(render);
  }

  return(r_val);
}

