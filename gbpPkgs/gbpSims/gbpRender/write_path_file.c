#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gd.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpRender.h>

void write_path_file(render_info *render,int frame){

  if(SID.I_am_Master){
     // Write header
     if(frame==0){
        char  filename[MAX_FILENAME_LENGTH];
        FILE *fp_out=NULL;

        // Write render file
        sprintf(filename,"%s/render.txt",render->filename_out_dir);
        fp_out=fopen(filename,"w");
        fprintf(fp_out,"%%n_interpolate          %d\n", render->n_interpolate);
        fprintf(fp_out,"%%n_frames               %d\n", render->n_frames);
        fprintf(fp_out,"%%filename_SSimPL_root   %s\n", render->filename_SSimPL_root);
        fprintf(fp_out,"%%filename_halos_version %s\n", render->filename_halos_version);
        fprintf(fp_out,"%%filename_trees_version %s\n", render->filename_trees_version);
        fprintf(fp_out,"%%snap_filename_root     %s\n", render->snap_filename_root);
        fprintf(fp_out,"%%mark_filename_root     %s\n", render->mark_filename_root);
        fprintf(fp_out,"%%smooth_filename_root   %s\n", render->smooth_filename_root);
        fprintf(fp_out,"%%snap_a_list_filename   %s\n", render->snap_a_list_filename);
        fprintf(fp_out,"%%flag_comoving          %d\n", render->flag_comoving);
        fprintf(fp_out,"%%flag_fade              %d\n", render->flag_fade);
        fprintf(fp_out,"%%alpha_fade             %le\n",render->alpha_fade);
        fprintf(fp_out,"%%flag_force_periodic    %d\n", render->flag_force_periodic);
        fprintf(fp_out,"%%flag_read_marked       %d\n", render->flag_read_marked);
        fprintf(fp_out,"%%flag_add_absorption    %d\n", render->flag_add_absorption);
        fprintf(fp_out,"%%f_absorption           %le\n",render->f_absorption);
        fprintf(fp_out,"%%w_mode                 %d\n", render->w_mode);
        fprintf(fp_out,"%%v_mode                 %d\n", render->v_mode);
        fclose(fp_out);

        // Write camera file
        sprintf(filename,"%s/camera.txt",render->filename_out_dir);
        fp_out=fopen(filename,"w");
        camera_info *camera=render->camera;
        fprintf(fp_out,"%%camera_mode         %d\n",     camera->camera_mode);
        fprintf(fp_out,"%%colour_table        %d\n",     camera->colour_table);
        fprintf(fp_out,"%%flag_velocity_space %d\n",     camera->flag_velocity_space);
        fprintf(fp_out,"%%width               %d\n",     camera->width);
        fprintf(fp_out,"%%height              %d\n",     camera->height);
        fprintf(fp_out,"%%stereo_ratio        %le\n",    camera->stereo_ratio);
        fprintf(fp_out,"%%f_near_field        %le\n",    camera->f_near_field);
        fprintf(fp_out,"%%f_taper_field       %le\n",    camera->f_taper_field);
        fprintf(fp_out,"%%f_image_plane       %le\n",    camera->f_image_plane);
        fprintf(fp_out,"%%RGB_mode            %d\n",     camera->RGB_mode);
        fprintf(fp_out,"%%flag_calc_Z_image   %d\n",     camera->flag_calc_Z_image);
        fprintf(fp_out,"%%RGB_param           %s\n",     camera->RGB_param);
        fprintf(fp_out,"%%RGB_range           %le %le\n",camera->RGB_range[0],camera->RGB_range[1]);
        fclose(fp_out);
     }

     // Write perspective entry
     char  filename[MAX_FILENAME_LENGTH];
     FILE *fp_out=NULL;
     sprintf(filename,"%s/perspective.txt",render->filename_out_dir);
     camera_info      *camera     =render->camera;
     perspective_info *perspective=camera->perspective;
     if(frame==0){
        int i_column=1;
        fp_out=fopen(filename,"w");
        fprintf(fp_out,"# Perspective file\n");
        fprintf(fp_out,"# Column (%02d): Frame\n",           i_column++);
        fprintf(fp_out,"#        (%02d): Snapshot\n",        i_column++);
        fprintf(fp_out,"#        (%02d): Expansion factor\n",i_column++);
        fprintf(fp_out,"#        (%02d): p_o_x\n",           i_column++);
        fprintf(fp_out,"#        (%02d): p_o_y\n",           i_column++);
        fprintf(fp_out,"#        (%02d): p_o_z\n",           i_column++);
        fprintf(fp_out,"#        (%02d): p_c_x\n",           i_column++);
        fprintf(fp_out,"#        (%02d): p_c_y\n",           i_column++);
        fprintf(fp_out,"#        (%02d): p_c_z\n",           i_column++);
        fprintf(fp_out,"#        (%02d): d_o\n",             i_column++);
        fprintf(fp_out,"#        (%02d): theta\n",           i_column++);
        fprintf(fp_out,"#        (%02d): phi\n",             i_column++);
        fprintf(fp_out,"#        (%02d): FOV\n",             i_column++);
        fprintf(fp_out,"#        (%02d): focus_shift_x\n",   i_column++);
        fprintf(fp_out,"#        (%02d): focus_shift_y\n",   i_column++);
     }
     else
        fp_out=fopen(filename,"a");
     int snap_best;
     pick_best_snap(perspective->time,render->snap_a_list,render->n_snap_a_list,&snap_best,NULL);
     fprintf(fp_out,"%04d %03d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n",
             frame,
             snap_best,
             perspective->time,
             perspective->p_o[0],
             perspective->p_o[1],
             perspective->p_o[2],
             perspective->p_c[0],
             perspective->p_c[1],
             perspective->p_c[2],
             perspective->d_o,
             perspective->theta,
             perspective->phi,
             perspective->FOV,
             perspective->focus_shift_x,
             perspective->focus_shift_y);
     fclose(fp_out);
  } 
}

