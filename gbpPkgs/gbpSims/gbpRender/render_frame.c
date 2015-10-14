#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>
#include <gbpRender.h>

void rotate_particle(double   x_hat,
                     double   y_hat,
                     double   z_hat,
                     double   theta,
                     GBPREAL *x_i,
                     GBPREAL *y_i,
                     GBPREAL *z_i);
void rotate_particle(double   x_hat,
                     double   y_hat,
                     double   z_hat,
                     double   theta,
                     GBPREAL *x_i,
                     GBPREAL *y_i,
                     GBPREAL *z_i){
  double X;
  double Y;
  double Z;
  double W;
  double x_p;
  double y_p;
  double z_p;
  X     =x_hat*sin(theta/2);
  Y     =y_hat*sin(theta/2);
  Z     =z_hat*sin(theta/2);
  W     =      cos(theta/2);
  x_p   =(double)(*x_i);
  y_p   =(double)(*y_i);
  z_p   =(double)(*z_i);
  (*x_i)=(1.-2.*Y*Y-2.*Z*Z)*x_p+(   2.*X*Y-2.*Z*W)*y_p+(   2.*X*Z+2.*Y*W)*z_p;
  (*y_i)=(   2.*X*Y+2.*Z*W)*x_p+(1.-2.*X*X-2.*Z*Z)*y_p+(   2.*Y*Z-2.*X*W)*z_p;
  (*z_i)=(   2.*X*Z-2.*Y*W)*x_p+(   2.*Y*Z+2.*X*W)*y_p+(1.-2.*X*X-2.*Y*Y)*z_p;
}

int set_pixel_space(float   h_i,
                    float   x_i,
                    float   y_i,
                    float   f_i,
                    double  xmin,
                    double  ymin,
                    double  FOV_x,
                    double  FOV_y,
                    double  pixel_size_x,
                    double  pixel_size_y,
                    double  radius_kernel_max,
                    double *radius2_norm,
                    double *radius_kernel,
                    double *part_pos_x,
                    double *part_pos_y,
                    int    *kx_min,
                    int    *kx_max,
                    int    *ky_min,
                    int    *ky_max);
int set_pixel_space(float   h_i,
                    float   x_i,
                    float   y_i,
                    float   f_i,
                    double  xmin,
                    double  ymin,
                    double  FOV_x,
                    double  FOV_y,
                    double  pixel_size_x,
                    double  pixel_size_y,
                    double  radius_kernel_max,
                    double *radius2_norm,
                    double *radius_kernel,
                    double *part_pos_x,
                    double *part_pos_y,
                    int    *kx_min,
                    int    *kx_max,
                    int    *ky_min,
                    int    *ky_max){
   double part_h_xy;
   part_h_xy       =(double)h_i*(double)f_i;
   (*part_pos_x)   =(double)x_i*(double)f_i;
   (*part_pos_y)   =(double)y_i*(double)f_i;
   (*radius2_norm) =1./(part_h_xy*part_h_xy);
   (*radius_kernel)=radius_kernel_max*part_h_xy;
   (*kx_min)       =(int)(((*part_pos_x)-(*radius_kernel)-xmin)/pixel_size_x);
   (*kx_max)       =(int)(((*part_pos_x)+(*radius_kernel)-xmin)/pixel_size_x+ONE_HALF);
   (*ky_min)       =(int)(((*part_pos_y)-(*radius_kernel)-ymin)/pixel_size_y);
   (*ky_max)       =(int)(((*part_pos_y)+(*radius_kernel)-ymin)/pixel_size_y+ONE_HALF);
}

void transform_particle(GBPREAL *x_i,
                        GBPREAL *y_i,
                        GBPREAL *z_i,
                        double   x_o,
                        double   y_o,
                        double   z_o,
                        double   x_hat,
                        double   y_hat,
                        double   z_hat,
                        double   d_o,
                        double   stereo_offset,
                        double   theta,
                        double   theta_roll,
                        double   box_size,
                        double   expansion_factor,
                        double   focus_shift_x,
                        double   focus_shift_y,
                        int      flag_comoving,
                        int      flag_force_periodic);
void transform_particle(GBPREAL *x_i,
                        GBPREAL *y_i,
                        GBPREAL *z_i,
                        double   x_o,
                        double   y_o,
                        double   z_o,
                        double   x_hat,
                        double   y_hat,
                        double   z_hat,
                        double   d_o,
                        double   stereo_offset,
                        double   theta,
                        double   theta_roll,
                        double   box_size,
                        double   expansion_factor,
                        double   focus_shift_x,
                        double   focus_shift_y,
                        int      flag_comoving,
                        int      flag_force_periodic){

   // Shift to object-centred coordinates
   (*x_i)-=(GBPREAL)x_o;
   (*y_i)-=(GBPREAL)y_o;
   (*z_i)-=(GBPREAL)z_o;

   // Centre the periodic box   
   if(!flag_comoving || flag_force_periodic){
     force_periodic(x_i,-0.5*(GBPREAL)box_size,(GBPREAL)box_size);
     force_periodic(y_i,-0.5*(GBPREAL)box_size,(GBPREAL)box_size);
     force_periodic(z_i,-0.5*(GBPREAL)box_size,(GBPREAL)box_size);
   }
   
   // Rotate about the origen to place camera position on z-axis
   rotate_particle(x_hat,
                   y_hat,
                   z_hat,
                   theta,
                   x_i,
                   y_i,
                   z_i);

   // Apply roll angle
   rotate_particle(0.,
                   0.,
                   1.,
                   (-theta_roll),
                   x_i,
                   y_i,
                   z_i);
   
   // Convert to proper coordinates?
   if(!flag_comoving){
     (*x_i)*=expansion_factor;
     (*y_i)*=expansion_factor;
     (*z_i)*=expansion_factor;
   }
   
   // Apply focus shift
   (*x_i)-=focus_shift_x;
   (*y_i)-=focus_shift_y;

   // Shift zero to the camera position
   (*z_i)+=d_o;

}

typedef struct map_quantities_info map_quantities_info;
struct map_quantities_info{
   int           flag_weigh;
   int           flag_line_integral;
   int          *ptype_used;
   size_t        n_particles;
   GBPREAL     **x;
   GBPREAL     **y;
   GBPREAL     **z;
   float       **h_smooth;
   float       **rho;
   float       **sigma;
   double        mass_array;
   interp_info  *transfer_rho;
   int           flag_transfer_rho_log;
   interp_info  *transfer_sigma;
   int           flag_transfer_sigma_log;
   int           flag_comoving;
   double        inv_expansion_factor_cubed;
   int           v_mode;
   int           w_mode;
};

void free_particle_map_quantities(map_quantities_info *mq);
void free_particle_map_quantities(map_quantities_info *mq){
   SID_free(SID_FARG mq->h_smooth);
   SID_free(SID_FARG mq->x);
   SID_free(SID_FARG mq->y);
   SID_free(SID_FARG mq->z);
   SID_free(SID_FARG mq->rho);
   SID_free(SID_FARG mq->sigma);
   SID_free(SID_FARG mq->ptype_used);
}

void init_particle_map_quantities(map_quantities_info *mq,render_info *render,ADaPS *transfer_list,int flag_comoving,double expansion_factor);
void init_particle_map_quantities(map_quantities_info *mq,render_info *render,ADaPS *transfer_list,int flag_comoving,double expansion_factor){
  int i_snap;

  SID_log("Initializing quantities...",SID_LOG_OPEN);

  // Defaults
  mq->flag_weigh             =FALSE;
  mq->flag_line_integral     =FALSE;
  mq->flag_transfer_sigma_log=FALSE;
  mq->flag_transfer_rho_log  =FALSE;
  mq->n_particles            =0;
  mq->mass_array             =0.;
  mq->h_smooth               =(float   **)SID_malloc(sizeof(float *)*render->n_interpolate);
  mq->x                      =(GBPREAL **)SID_malloc(sizeof(float *)*render->n_interpolate);
  mq->y                      =(GBPREAL **)SID_malloc(sizeof(float *)*render->n_interpolate);
  mq->z                      =(GBPREAL **)SID_malloc(sizeof(float *)*render->n_interpolate);
  mq->rho                    =(float   **)SID_malloc(sizeof(float *)*render->n_interpolate);
  mq->sigma                  =(float   **)SID_malloc(sizeof(float *)*render->n_interpolate);
  mq->transfer_sigma         =NULL;
  mq->transfer_rho           =NULL;

  // Comoving?
  mq->flag_comoving=flag_comoving;
  if(mq->flag_comoving)
     mq->inv_expansion_factor_cubed=1/(expansion_factor*expansion_factor*expansion_factor);
  else
     mq->inv_expansion_factor_cubed=1.;

  // Initialize the array of used particle types
  int i_type;
  mq->ptype_used=(int *)SID_malloc(sizeof(int)*N_GADGET_TYPE);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
     mq->ptype_used[i_type]=FALSE;

  // Render a dark matter velocity dispersion map
  if(!strcmp(render->camera->RGB_param,"sigma_v_dark") && !strcmp(render->camera->Y_param,"tau_dark")){
    mq->ptype_used[GADGET_TYPE_DARK]=TRUE;
    mq->flag_weigh                  =TRUE;
    mq->flag_line_integral          =TRUE;
    mq->n_particles                 =((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_dark"))[0];
    mq->v_mode=MAKE_MAP_MODE_SIGMA;
    mq->w_mode=MAKE_MAP_MODE_RHO;
    if(ADaPS_exist(render->plist_list[0]->data,"rho_dark")){
      for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
         mq->rho[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"rho_dark");
      if(ADaPS_exist(transfer_list,"rho_dark")){
        if(ADaPS_exist(transfer_list,"rho_dark_log")) 
          mq->flag_transfer_rho_log=TRUE;
        mq->transfer_rho=(interp_info *)ADaPS_fetch(transfer_list,"rho_dark");
      }
      else{
        mq->flag_transfer_rho_log=FALSE;
        mq->transfer_rho         =NULL;
      }
    }
    else
      SID_trap_error("No densities available in make_map.",ERROR_LOGIC);

    // Use sigma_v for values
    if(ADaPS_exist(render->plist_list[0]->data,"sigma_v_dark")){
      for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
         mq->sigma[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"sigma_v_dark");
      if(ADaPS_exist(transfer_list,"sigma_v_dark")){
        if(ADaPS_exist(transfer_list,"sigma_v_dark_log")) 
          mq->flag_transfer_sigma_log=TRUE;
        mq->transfer_sigma=(interp_info *)ADaPS_fetch(transfer_list,"sigma_v_dark");
      } 
      else{
        mq->flag_transfer_sigma_log=FALSE;
        mq->transfer_sigma         =NULL;
      }
    }
    else
      SID_trap_error("No sigma_v's available in make_map.",ERROR_LOGIC);
  }
  else if(!strcmp(render->camera->Y_param,"sigma_v_dark") && !strcmp(render->camera->RGB_param,"tau_dark")){
    mq->ptype_used[GADGET_TYPE_DARK]=TRUE;
    mq->flag_weigh                  =TRUE;
    mq->flag_line_integral          =TRUE;
    mq->n_particles                 =((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_dark"))[0];
    mq->v_mode=MAKE_MAP_MODE_SIGMA;
    mq->w_mode=MAKE_MAP_MODE_RHO;
    if(ADaPS_exist(render->plist_list[0]->data,"sigma_v_dark")){
      for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
         mq->rho[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"sigma_v_dark");
      if(ADaPS_exist(transfer_list,"sigma_v_dark")){
        if(ADaPS_exist(transfer_list,"sigma_v_dark_log"))
          mq->flag_transfer_rho_log=TRUE;
        mq->transfer_rho=(interp_info *)ADaPS_fetch(transfer_list,"sigma_v_dark");
      }
      else{
        mq->flag_transfer_rho_log=FALSE;
        mq->transfer_rho         =NULL;
      }
    }
    else
      SID_trap_error("No densities available in make_map.",ERROR_LOGIC);

    // Use rho for values
    if(ADaPS_exist(render->plist_list[0]->data,"rho_dark")){
      for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
         mq->sigma[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"rho_dark");
      if(ADaPS_exist(transfer_list,"rho_dark")){
        if(ADaPS_exist(transfer_list,"rho_dark_log"))
          mq->flag_transfer_sigma_log=TRUE;
        mq->transfer_sigma=(interp_info *)ADaPS_fetch(transfer_list,"rho_dark");
      }
      else{
        mq->flag_transfer_sigma_log=FALSE;
        mq->transfer_sigma         =NULL;
      }
    }
    else
      SID_trap_error("No sigma_v's available in make_map.",ERROR_LOGIC);
  }
  else if(!strcmp(render->camera->RGB_param,"sigma_v_dark") && render->camera->flag_velocity_space){
    mq->ptype_used[GADGET_TYPE_DARK]=TRUE;
    mq->flag_weigh                  =FALSE;
    mq->flag_line_integral          =TRUE;
    mq->n_particles                 =((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_dark"))[0];
    mq->v_mode=MAKE_MAP_MODE_SIGMA;
    mq->w_mode=MAKE_MAP_INV_SIGMA;

    // Use sigma_v for values
    if(ADaPS_exist(render->plist_list[0]->data,"sigma_v_dark")){
      for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
         mq->sigma[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"sigma_v_dark");
      if(ADaPS_exist(transfer_list,"sigma_v_dark")){
        if(ADaPS_exist(transfer_list,"sigma_v_dark_log"))
          mq->flag_transfer_sigma_log=TRUE;
        mq->transfer_sigma=(interp_info *)ADaPS_fetch(transfer_list,"sigma_v_dark");
      }
      else{
        mq->flag_transfer_sigma_log=FALSE;
        mq->transfer_sigma         =NULL;
      }
    }
    else
      SID_trap_error("No sigma_v's available in make_map.",ERROR_LOGIC);
  }
  else
    SID_trap_error("Unknown rendering configuration RGB={%s} Y={%s} w/ flag_velocity_space=%d.",ERROR_LOGIC,
                   render->camera->RGB_param,
                   render->camera->Y_param,
                   render->camera->flag_velocity_space);

  // Initialize arrays
  int n_type_used=0;
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
     if(mq->ptype_used[i_type]){
        // This code can use one-and-only-one particle type at a time at the moment
        n_type_used++;
        if(n_type_used>1)
           SID_trap_error("An invalid number of particle types (%d) are being used in make_map.",ERROR_LOGIC,n_type_used);
        if(render->camera->flag_velocity_space){
           for(i_snap=0;i_snap<render->n_interpolate;i_snap++){
              mq->x[i_snap]=(GBPREAL *)ADaPS_fetch(render->plist_list[i_snap]->data,"vx_%s",render->plist_list[i_snap]->species[i_type]);
              mq->y[i_snap]=(GBPREAL *)ADaPS_fetch(render->plist_list[i_snap]->data,"vy_%s",render->plist_list[i_snap]->species[i_type]);
              mq->z[i_snap]=(GBPREAL *)ADaPS_fetch(render->plist_list[i_snap]->data,"vz_%s",render->plist_list[i_snap]->species[i_type]);
           }
           for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
              mq->h_smooth[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"sigma_v_%s",render->plist_list[i_snap]->species[i_type]);
        }
        else{
           for(i_snap=0;i_snap<render->n_interpolate;i_snap++){
              mq->x[i_snap]=(GBPREAL *)ADaPS_fetch(render->plist_list[i_snap]->data,"x_%s",render->plist_list[i_snap]->species[i_type]);
              mq->y[i_snap]=(GBPREAL *)ADaPS_fetch(render->plist_list[i_snap]->data,"y_%s",render->plist_list[i_snap]->species[i_type]);
              mq->z[i_snap]=(GBPREAL *)ADaPS_fetch(render->plist_list[i_snap]->data,"z_%s",render->plist_list[i_snap]->species[i_type]);
           }
           for(i_snap=0;i_snap<render->n_interpolate;i_snap++)
              mq->h_smooth[i_snap]=(float *)ADaPS_fetch(render->plist_list[i_snap]->data,"r_smooth_%s",render->plist_list[i_snap]->species[i_type]);
        }
     }
  }
  SID_log("%zd eligible for rendering...Done.",SID_LOG_CLOSE,mq->n_particles);

}

void set_particle_map_quantities(render_info *render,map_quantities_info *mq,int mode,size_t i_particle,
                                 float    box_size,
                                 float    half_box,
                                 float   *x_i,
                                 float   *y_i,
                                 float   *z_i,
                                 float   *h_i,
                                 float   *v_i,
                                 float   *w_i);
void set_particle_map_quantities(render_info *render,map_quantities_info *mq,int mode,size_t i_particle,
                                 float   box_size,
                                 float   half_box,
                                 float  *x_i,
                                 float  *y_i,
                                 float  *z_i,
                                 float  *h_i,
                                 float  *v_i,
                                 float  *w_i){
  if(render->n_interpolate==1){ 
     (*x_i)=mq->x[0][i_particle];
     (*y_i)=mq->y[0][i_particle];
     (*z_i)=mq->z[0][i_particle];
  }
  else if(render->n_interpolate==2){
     float dx=(float)(mq->x[1][i_particle]-mq->x[0][i_particle]);
     float dy=(float)(mq->y[1][i_particle]-mq->y[0][i_particle]);
     float dz=(float)(mq->z[1][i_particle]-mq->z[0][i_particle]);
     if(dx> half_box) dx-=box_size;
     if(dx<-half_box) dx+=box_size;
     if(dy> half_box) dy-=box_size;
     if(dy<-half_box) dy+=box_size;
     if(dz> half_box) dz-=box_size;
     if(dz<-half_box) dz+=box_size;
     (*x_i)=mq->x[0][i_particle]+render->f_interpolate*dx;
     (*y_i)=mq->y[0][i_particle]+render->f_interpolate*dy;
     (*z_i)=mq->z[0][i_particle]+render->f_interpolate*dz;
  }
  else
    SID_trap_error("n_interpolate>2 not supported (yet).",ERROR_LOGIC);
  if(mode==TRUE){
     if(render->n_interpolate==1)
        (*h_i)=mq->h_smooth[0][i_particle];
     else if(render->n_interpolate==2)
        (*h_i)=mq->h_smooth[0][i_particle]+render->f_interpolate*(mq->h_smooth[1][i_particle]-mq->h_smooth[0][i_particle]);
     // Set the particle weighting
     switch(mq->w_mode){
        case MAKE_MAP_MODE_RHO:
           if(render->n_interpolate==1)
              (*w_i)=mq->rho[0][i_particle];
           else if(render->n_interpolate==2){
              double rho_0=mq->rho[0][i_particle];
              double rho_1=mq->rho[1][i_particle];
              (*w_i)=take_alog10(rho_0+render->f_interpolate*(rho_1-rho_0));
           }
           if(mq->transfer_rho!=NULL){
              switch(mq->flag_transfer_rho_log){
                 case TRUE:
                    (*w_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,take_log10((double)(*w_i)))));
                    break;
                 case FALSE:
                    (*w_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,(double)(*w_i))));
                    break;
              }
           }
           /*
           if(!mq->flag_comoving)
              (*w_i)*=mq->inv_expansion_factor_cubed;
           */
           break;
        case MAKE_MAP_INV_SIGMA:
           if(render->n_interpolate==1)
              (*w_i)=1./mq->sigma[0][i_particle];
           else if(render->n_interpolate==2){
              double sigma_0=1./mq->sigma[0][i_particle];
              double sigma_1=1./mq->sigma[1][i_particle];
              (*w_i)=take_alog10(sigma_0+render->f_interpolate*(sigma_1-sigma_0));
           }
           if(mq->transfer_sigma!=NULL){
              switch(mq->flag_transfer_sigma_log){
                 case TRUE:
                    (*w_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_sigma,take_log10((double)(*w_i)))));
                    break;
                 case FALSE:
                    (*w_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_sigma,(double)(*w_i))));
                    break;
              }
           }
           break;
        case MAKE_MAP_NO_WEIGHTING:
           (*w_i)=1.;
           break;
        default:
           SID_trap_error("Unknown w_mode (%d) in make_map.",ERROR_LOGIC,mq->v_mode);
           break;
     }
     // Set the particle value
     switch(mq->v_mode){
        case MAKE_MAP_MODE_SIGMA:
           if(render->n_interpolate==1)
              (*v_i)=mq->sigma[0][i_particle];
           else if(render->n_interpolate==2){
              double sigma_0=mq->sigma[0][i_particle];
              double sigma_1=mq->sigma[1][i_particle];
              (*v_i)=take_alog10(sigma_0+render->f_interpolate*(sigma_1-sigma_0));
           }
           if(mq->transfer_sigma!=NULL){
              switch(mq->flag_transfer_sigma_log){
                 case TRUE:
                    (*v_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,take_log10((double)(*v_i)))));
                    break;
                 case FALSE:
                    (*v_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,(double)(*v_i))));
                    break;
              }
           }
           break;
        default:
           SID_trap_error("Unknown v_mode (%d) in make_map.",ERROR_LOGIC,mq->v_mode);
           break;
     }
  }
}

float compute_f_stretch(double d_image_plane,float z_i,int flag_plane_parallel);
float compute_f_stretch(double d_image_plane,float z_i,int flag_plane_parallel){
   float f_i;
   switch(flag_plane_parallel){
      case FALSE:
         if(z_i>0.)
            f_i=(float)d_image_plane/(float)z_i;
         else
            f_i=0.;
         break;
      case TRUE:
         f_i=1.;
         break;
   }
   return(f_i);
}

// Compute the angle and the axis of rotation
//   needed to place the camera at (0,0,-d_o)
//   with the object at (0,0,0)
void compute_perspective_transformation(double  x_o,
                                        double  y_o,
                                        double  z_o,
                                        double  x_c,
                                        double  y_c,
                                        double  z_c,
                                        double  unit_factor,
                                        const char *unit_text,
                                        double  f_image_plane,
                                        double  stereo_offset,
                                        double *FOV_x,
                                        double *FOV_y,
                                        double *d_o,
                                        double *x_o_out,
                                        double *y_o_out,
                                        double *z_o_out,
                                        double *x_c_out,
                                        double *y_c_out,
                                        double *z_c_out,
                                        double *x_hat,
                                        double *y_hat,
                                        double *z_hat,
                                        double *theta,
                                        double *theta_roll);
void compute_perspective_transformation(double  x_o,
                                        double  y_o,
                                        double  z_o,
                                        double  x_c,
                                        double  y_c,
                                        double  z_c,
                                        double  unit_factor,
                                        const char *unit_text,
                                        double  f_image_plane,
                                        double  stereo_offset,
                                        double *FOV_x,
                                        double *FOV_y,
                                        double *d_o,
                                        double *x_o_out,
                                        double *y_o_out,
                                        double *z_o_out,
                                        double *x_c_out,
                                        double *y_c_out,
                                        double *z_c_out,
                                        double *x_hat,
                                        double *y_hat,
                                        double *z_hat,
                                        double *theta,
                                        double *theta_roll){
  double d_x_o;
  double d_y_o;
  double d_z_o;
  double d_xy;
  d_x_o   = x_o-x_c;
  d_y_o   = y_o-y_c;
  d_z_o   = z_o-z_c;
  d_xy    = sqrt(pow(d_x_o,2.)+pow(d_y_o,2.));
  (*d_o)  = sqrt(pow(d_x_o,2.)+pow(d_y_o,2.)+pow(d_z_o,2.));
  (*x_hat)= d_y_o/d_xy;
  (*y_hat)=-d_x_o/d_xy;
  (*z_hat)= 0.;
  (*theta)= acos(d_z_o/(*d_o));
  if(sqrt(d_x_o*d_x_o+d_y_o*d_y_o)>0.){
    (*theta_roll)=acos(-d_y_o/sqrt(d_x_o*d_x_o+d_y_o*d_y_o));
    if(d_x_o<0.)
      (*theta_roll)=TWO_PI-(*theta_roll);
  }
  else
    (*theta_roll)=0.;

  // Apply stereo offsets
  (*x_o_out)=x_o;
  (*y_o_out)=y_o;
  (*z_o_out)=z_o;
  (*x_c_out)=x_c;
  (*y_c_out)=y_c;
  (*z_c_out)=z_c;
  if(stereo_offset!=0.){
     double Dx_stereo;
     double Dy_stereo;
     double Dz_stereo;
     double theta_roll_stereo=0.;
     double cos_theta_over_d_xy;
     SID_log("Forcing theta=0 in stereo projection.",SID_LOG_COMMENT);
     if(d_x_o>0.){
        if(d_y_o>0.){
           Dx_stereo=-stereo_offset*cos(theta_roll_stereo)*fabs(d_y_o)/d_xy;
           Dy_stereo= stereo_offset*cos(theta_roll_stereo)*fabs(d_x_o)/d_xy;
        }
        else{
           Dx_stereo= stereo_offset*cos(theta_roll_stereo)*fabs(d_y_o)/d_xy;
           Dy_stereo= stereo_offset*cos(theta_roll_stereo)*fabs(d_x_o)/d_xy; 
        }
     }
     else{
        if(d_y_o>0.){
           Dx_stereo=-stereo_offset*cos(theta_roll_stereo)*fabs(d_y_o)/d_xy; 
           Dy_stereo=-stereo_offset*cos(theta_roll_stereo)*fabs(d_x_o)/d_xy; 
        }
        else{
           Dx_stereo= stereo_offset*cos(theta_roll_stereo)*fabs(d_y_o)/d_xy; 
           Dy_stereo=-stereo_offset*cos(theta_roll_stereo)*fabs(d_x_o)/d_xy; 
        }
     }
     if(theta_roll_stereo>0.)
        Dz_stereo=sqrt(stereo_offset*stereo_offset-Dx_stereo*Dx_stereo-Dy_stereo*Dy_stereo);
     else if(theta_roll_stereo==0.)
        Dz_stereo=0.;
     else
        Dz_stereo=-sqrt(stereo_offset*stereo_offset-Dx_stereo*Dx_stereo-Dy_stereo*Dy_stereo);
     if(stereo_offset!=0){
        SID_log("Stereo offset results: (offset=%10.3le [%s])",SID_LOG_OPEN,   stereo_offset/unit_factor,unit_text);
        SID_log("(x_o,x_c,D_x)=(%10.3le,%10.3le,%10.3le) [%s]",SID_LOG_COMMENT,x_o/unit_factor,x_c/unit_factor,Dx_stereo/unit_factor,unit_text);
        SID_log("(y_o,y_c,D_y)=(%10.3le,%10.3le,%10.3le) [%s]",SID_LOG_COMMENT,y_o/unit_factor,y_c/unit_factor,Dy_stereo/unit_factor,unit_text);
        SID_log("(z_o,z_c,D_z)=(%10.3le,%10.3le,%10.3le) [%s]",SID_LOG_COMMENT,z_o/unit_factor,z_c/unit_factor,Dz_stereo/unit_factor,unit_text);
        SID_log("",SID_LOG_SILENT_CLOSE);
     }
     (*x_o_out)+=Dx_stereo;
     (*y_o_out)+=Dy_stereo;
     (*z_o_out)+=Dz_stereo;
     (*x_c_out)+=Dx_stereo;
     (*y_c_out)+=Dy_stereo;
     (*z_c_out)+=Dz_stereo;
  }

  // Recompute d_o to place the object in front/behind of the image plane
  (*d_o)  /=f_image_plane;
  (*FOV_x)*=f_image_plane;
  (*FOV_y)*=f_image_plane;

}

int check_if_particle_marked(char **mark,int i_species,size_t i_particle,char *c_i);
int check_if_particle_marked(char **mark,int i_species,size_t i_particle,char *c_i){
   // Return FALSE only if mark[i_species][i_particle] is defined ...
   if(mark[i_species]==NULL){
      (*c_i)=1; // default to white
      return(TRUE);
   }
   // ... and equal to FALSE.
   else{
      (*c_i)=mark[i_species][i_particle];
      return(!(mark[i_species][i_particle]==FALSE));
   }
}

void init_make_map_noabs(render_info *render,
                         double       x_o,
                         double       y_o,
                         double       z_o,
                         double       x_c,
                         double       y_c,
                         double       z_c,
                         double       unit_factor,
                         const char  *unit_text,
                         double       f_image_plane,
                         double       box_size,
                         double       FOV_x_in,
                         double       FOV_y_in,
                         double       xmin,
                         double       ymin,
                         double       pixel_size_x,
                         double       pixel_size_y,
                         double       radius_kernel_max,
                         int          nx,
                         int          ny,
                         double       expansion_factor,
                         double       focus_shift_x,
                         double       focus_shift_y,
                         double       d_near_field,
                         double       stereo_offset,
                         int          flag_comoving,
                         int          flag_force_periodic,
                         int          camera_mode,
                         int         *flag_weigh,
                         int         *flag_line_integral,
                         float       **x,
                         float       **y,
                         float       **z,
                         float       **h_smooth,
                         float       **f_stretch,
                         float       **value,
                         float       **weight,
                         char        **colour,
                         size_t      **z_index,
                         int          *i_x_min_local_return,
                         int          *i_x_max_local_return,
                         size_t       *n_particles);
void init_make_map_noabs(render_info *render,
                         double       x_o,
                         double       y_o,
                         double       z_o,
                         double       x_c,
                         double       y_c,
                         double       z_c,
                         double       unit_factor,
                         const char  *unit_text,
                         double       f_image_plane,
                         double       box_size,
                         double       FOV_x_in,
                         double       FOV_y_in,
                         double       xmin,
                         double       ymin,
                         double       pixel_size_x,
                         double       pixel_size_y,
                         double       radius_kernel_max,
                         int          nx,
                         int          ny,
                         double       expansion_factor,
                         double       focus_shift_x,
                         double       focus_shift_y,
                         double       d_near_field,
                         double       stereo_offset,
                         int          flag_comoving,
                         int          flag_force_periodic,
                         int          camera_mode,
                         int         *flag_weigh,
                         int         *flag_line_integral,
                         float       **x,
                         float       **y,
                         float       **z,
                         float       **h_smooth,
                         float       **f_stretch,
                         float       **value,
                         float       **weight,
                         char        **colour,
                         size_t      **z_index,
                         int          *i_x_min_local_return,
                         int          *i_x_max_local_return,
                         size_t       *n_particles){
  int     *ptype_used;
  int      i_type;
  float   *rho;
  double  *mass;
  float   *sigma_v;
  size_t   i_particle;
  size_t   j_particle;
  size_t   k_particle;
  double   mass_array;
  size_t   n_particles_species;
  float   *x_temp;
  float   *y_temp;
  float   *z_temp;
  float   *h_smooth_temp;
  double   x_hat;
  double   y_hat;
  double   z_hat;
  double   theta;
  double   theta_roll;
  double   d_x_o;
  double   d_y_o;
  double   d_z_o;
  double   d_o;
  double   d_hat;
  double   particle_radius;
  double   x_tmp,y_tmp,z_tmp;
  int      flag_log;
  double   z_test;
  int      flag_use_Gadget;
  float    box_size_float=(float)box_size;
  float    half_box=0.5*box_size_float;
  
  SID_log("Initializing projection-space...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Plane parallel projection?
  int flag_plane_parallel;
  if(check_mode_for_flag(camera_mode,CAMERA_PLANE_PARALLEL))
     flag_plane_parallel=TRUE;
  else
     flag_plane_parallel=FALSE;

  // Report the state of some flags
  if(flag_plane_parallel) SID_log("Plane-parallel  projection: ON", SID_LOG_COMMENT);
  else                    SID_log("Plane-parallel  projection: OFF",SID_LOG_COMMENT);
  if(flag_comoving)       SID_log("Comoving        projection: ON", SID_LOG_COMMENT);
  else                    SID_log("Comoving        projection: OFF",SID_LOG_COMMENT);
  if(flag_force_periodic) SID_log("Centre periodic projection: ON", SID_LOG_COMMENT);
  else                    SID_log("Centre periodic projection: OFF",SID_LOG_COMMENT);

  // Initialize the mapping quantities
  map_quantities_info mq;
  init_particle_map_quantities(&mq,render,render->camera->transfer_list,flag_comoving,expansion_factor);
  (*n_particles)       =mq.n_particles;
  (*flag_weigh)        =mq.flag_weigh;
  (*flag_line_integral)=mq.flag_line_integral;
  ptype_used           =mq.ptype_used;

  double x_o_in=x_o;
  double y_o_in=y_o;
  double z_o_in=z_o;
  double x_c_in=x_c;
  double y_c_in=y_c;
  double z_c_in=z_c;
  double FOV_x =FOV_x_in;
  double FOV_y =FOV_y_in;
  compute_perspective_transformation(x_o_in,
                                     y_o_in,
                                     z_o_in,
                                     x_c_in,
                                     y_c_in,
                                     z_c_in,
                                     unit_factor,
                                     unit_text,
                                     f_image_plane,
                                     stereo_offset,
                                     &FOV_x,
                                     &FOV_y,
                                     &d_o,
                                     &x_o,
                                     &y_o,
                                     &z_o,
                                     &x_c,
                                     &y_c,
                                     &z_c,
                                     &x_hat,
                                     &y_hat,
                                     &z_hat,
                                     &theta,
                                     &theta_roll);

  // The previous line sets d_o to the object distance.  Compute the distance to the image plane.
  double d_image_plane=d_o*f_image_plane;

  // Set mark arrays
  int    flag_mark_on=FALSE;
  char **mark        =(char **)SID_malloc(sizeof(char *)*N_GADGET_TYPE);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(ADaPS_exist(render->plist_list[0]->data,"mark_%s",render->plist_list[0]->species[i_type])){
         mark[i_type]=(char *)ADaPS_fetch(render->plist_list[0]->data,"mark_%s",render->plist_list[0]->species[i_type]);
         flag_mark_on=TRUE;
      }
      else
         mark[i_type]=NULL;
  }

  // Determine how many particles are contributing to each 
  //    column of the image
  float   x_i;
  float   y_i;
  float   z_i;
  float   h_i;
  float   f_i;
  float   v_i;
  float   w_i;
  char    c_i;
  int     i_x;

  // Set the range of columns assigned to each rank (all in this case)
  int    i_x_min_local=0;
  int    i_x_max_local=nx-1;

  // Count the number of particles contributing to each rank
  int     i_rank;
  size_t *n_rank;
  size_t *n_rank_local;
  size_t  n_particles_local;
  size_t  n_particles_visible_local=0;
  size_t  l_particle;
  n_rank      =(size_t *)SID_calloc(sizeof(size_t)*SID.n_proc);
  n_rank_local=(size_t *)SID_calloc(sizeof(size_t)*SID.n_proc);
  SID_log("Count particles...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
     for(i_type=0,j_particle=0,l_particle=0;i_type<N_GADGET_TYPE;i_type++){
       if(ptype_used[i_type] && ADaPS_exist(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type])){
         n_particles_species=((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type]))[0];
         for(i_particle=i_rank,k_particle=j_particle+i_rank;i_particle<n_particles_species;i_particle+=SID.n_proc,k_particle+=SID.n_proc){
            if(check_if_particle_marked(mark,i_type,i_particle,&c_i)){
   
               // Set the preoperties of the particle to be mapped (mode is FALSE because we DON'T need the angular size of the particle)
               set_particle_map_quantities(render,&mq,FALSE,k_particle,box_size_float,half_box,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);
   
               // Transform particle to render-coordinates
               transform_particle(&x_i,
                                  &y_i,
                                  &z_i,
                                  x_o,
                                  y_o,
                                  z_o,
                                  x_hat,
                                  y_hat,
                                  z_hat,
                                  d_o,
                                  stereo_offset,
                                  theta,
                                  theta_roll,
                                  box_size,
                                  expansion_factor,
                                  focus_shift_x,
                                  focus_shift_y,
                                  flag_comoving,
                                  flag_force_periodic);
   
               if(z_i>d_near_field){
                  n_rank_local[i_rank]++;
                  n_particles_visible_local++;
               }
            }
         }
         j_particle+=n_particles_species;
       }
     }
  }
  SID_Allreduce(n_rank_local,n_rank,SID.n_proc,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  n_particles_local=n_rank[SID.My_rank];
  (*n_particles)   =n_particles_local;
  SID_log("Done.",SID_LOG_CLOSE);

  // Report decomposition results
  if(SID.n_proc>1){
     SID_log("Results of domain decomposition:",SID_LOG_OPEN);
     for(i_rank=0;i_rank<SID.n_proc;i_rank++)
        SID_log("Rank %03d n_particles=%zd",SID_LOG_COMMENT,i_rank,n_rank[i_rank]);
     SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
  }

  // Determine the needed size of the comm buffers
  size_t n_buffer;
  size_t n_buffer_max;
  calc_max(n_rank_local,&n_buffer,SID.n_proc,SID_SIZE_T,CALC_MODE_DEFAULT);
  SID_Allreduce(&n_buffer,&n_buffer_max,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);

  // Broadcast particles
  size_t particle_index;
  SID_log("Broadcast particles...(buffer=%zd particles)...",SID_LOG_OPEN|SID_LOG_TIMER,n_buffer_max);
  float  *x_buffer;
  float  *y_buffer;
  float  *z_buffer;
  float  *h_buffer;
  float  *f_buffer;
  float  *v_buffer;
  float  *w_buffer;
  char   *c_buffer;
  (*x)        =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*y)        =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*z)        =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*h_smooth) =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*f_stretch)=(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*value)    =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*weight)   =(float *)SID_malloc(sizeof(float)*n_particles_local);
  x_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  y_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  z_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  h_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  f_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  v_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  w_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  if(flag_mark_on){
     (*colour)=(char  *)SID_malloc(sizeof(char)*n_particles_local);
     c_buffer =(char  *)SID_malloc(sizeof(char)*n_buffer);
  }
  else{
     (*colour)=NULL;
     c_buffer =NULL;
  }
  size_t i_particle_rank;
  size_t j_particle_rank=0;
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
     size_t particle_index=0;
     int    rank_to;
     int    rank_from;
     set_exchange_ring_ranks(&rank_to,&rank_from,i_rank);
     for(i_type=0,j_particle=0,l_particle=0;i_type<N_GADGET_TYPE;i_type++){
       if(ptype_used[i_type]){
         n_particles_species=((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type]))[0];
         for(i_particle=rank_to,k_particle=j_particle+rank_to;i_particle<n_particles_species;i_particle+=SID.n_proc,k_particle+=SID.n_proc){
            if(check_if_particle_marked(mark,i_type,i_particle,&c_i)){
               // Set the properties of the particle to be mapped 
               set_particle_map_quantities(render,&mq,TRUE,k_particle,box_size_float,half_box,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);
   
               // Transform particle to render-coordinates
               transform_particle(&x_i,
                                  &y_i,
                                  &z_i,
                                  x_o,
                                  y_o,
                                  z_o,
                                  x_hat,
                                  y_hat,
                                  z_hat,
                                  d_o,
                                  stereo_offset,
                                  theta,
                                  theta_roll,
                                  box_size,
                                  expansion_factor,
                                  focus_shift_x,
                                  focus_shift_y,
                                  flag_comoving,
                                  flag_force_periodic);
   
               if(z_i>d_near_field){
                  f_i=(float)compute_f_stretch(d_image_plane,z_i,flag_plane_parallel);
                  x_buffer[particle_index]=x_i;
                  y_buffer[particle_index]=y_i;
                  z_buffer[particle_index]=z_i;
                  h_buffer[particle_index]=h_i;
                  f_buffer[particle_index]=f_i;
                  w_buffer[particle_index]=w_i;
                  v_buffer[particle_index]=v_i;
                  if(c_buffer!=NULL)
                     c_buffer[particle_index]=(char)c_i;
                  particle_index++;
               }
            }
         }
         j_particle+=n_particles_species;
       }
     } // loop over particles
     if(particle_index!=n_rank_local[rank_to])
        SID_trap_error("Buffer particle count does not equal desired exchange count (ie. %zd!=%zd)",ERROR_LOGIC,particle_index,n_rank_local[rank_to]);

     // Perform exchanges
     size_t n_exchange;
     exchange_ring_buffer(x_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*x)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(y_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*y)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(z_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*z)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(h_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*h_smooth)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(f_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*f_stretch)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(w_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*weight)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(v_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*value)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     if(c_buffer!=NULL)
        exchange_ring_buffer(c_buffer,
                             sizeof(char),
                             n_rank_local[rank_to],
                             &((*colour)[j_particle_rank]),
                             &n_exchange,
                             i_rank);
     j_particle_rank+=n_exchange;
  } // i_rank

  // Sanity checks
  if(j_particle_rank!=n_particles_local)
     SID_trap_error("The wrong number of particles were received (ie. %zd!=%zd) on rank %d.",ERROR_LOGIC,
                    j_particle_rank,n_particles_local,SID.My_rank);

  SID_log("Done.",SID_LOG_CLOSE);

  // Sort the local particles by position
  merge_sort((*z),(size_t)(*n_particles),z_index,SID_FLOAT,SORT_COMPUTE_INDEX,FALSE);

  // Clean-up
  SID_free(SID_FARG n_rank);
  SID_free(SID_FARG n_rank_local);
  SID_free(SID_FARG x_buffer);
  SID_free(SID_FARG y_buffer);
  SID_free(SID_FARG z_buffer);
  SID_free(SID_FARG h_buffer);
  SID_free(SID_FARG f_buffer);
  SID_free(SID_FARG v_buffer);
  SID_free(SID_FARG w_buffer);
  SID_free(SID_FARG c_buffer);
  free_particle_map_quantities(&mq);
  SID_free(SID_FARG mark);

  (*i_x_min_local_return)=i_x_min_local;
  (*i_x_max_local_return)=i_x_max_local;

  SID_log("Done.",SID_LOG_CLOSE);
}

void init_make_map_abs(render_info *render,
                       double       x_o,
                       double       y_o,
                       double       z_o,
                       double       x_c,
                       double       y_c,
                       double       z_c,
                       double       unit_factor,
                       const char  *unit_text,
                       double       f_image_plane,
                       double       box_size,
                       double       FOV_x_in,
                       double       FOV_y_in,
                       double       xmin,
                       double       ymin,
                       double       pixel_size_x,
                       double       pixel_size_y,
                       double       radius_kernel_max,
                       int          nx,
                       int          ny,
                       double       expansion_factor,
                       double       focus_shift_x,
                       double       focus_shift_y,
                       double       d_near_field,
                       double       stereo_offset,
                       int          flag_comoving,
                       int          flag_force_periodic,
                       int          camera_mode,
                       int         *flag_weigh,
                       int         *flag_line_integral,
                       float       **x,
                       float       **y,
                       float       **z,
                       float       **h_smooth,
                       float       **f_stretch,
                       float       **value,
                       float       **weight,
                       char        **colour,
                       size_t      **z_index,
                       int          *i_x_min_local_return,
                       int          *i_x_max_local_return,
                       size_t       *n_particles);
void init_make_map_abs(render_info *render,
                       double       x_o,
                       double       y_o,
                       double       z_o,
                       double       x_c,
                       double       y_c,
                       double       z_c,
                       double       unit_factor,
                       const char  *unit_text,
                       double       f_image_plane,
                       double       box_size,
                       double       FOV_x_in,
                       double       FOV_y_in,
                       double       xmin,
                       double       ymin,
                       double       pixel_size_x,
                       double       pixel_size_y,
                       double       radius_kernel_max,
                       int          nx,
                       int          ny,
                       double       expansion_factor,
                       double       focus_shift_x,
                       double       focus_shift_y,
                       double       d_near_field,
                       double       stereo_offset,
                       int          flag_comoving,
                       int          flag_force_periodic,
                       int          camera_mode,
                       int         *flag_weigh,
                       int         *flag_line_integral,
                       float       **x,
                       float       **y,
                       float       **z,
                       float       **h_smooth,
                       float       **f_stretch,
                       float       **value,
                       float       **weight,
                       char        **colour,
                       size_t      **z_index,
                       int          *i_x_min_local_return,
                       int          *i_x_max_local_return,
                       size_t       *n_particles){
  int     *ptype_used;
  int      i_type;
  float   *rho;
  double  *mass;
  float   *sigma_v;
  size_t   i_particle;
  size_t   j_particle;
  size_t   k_particle;
  double   mass_array;
  size_t   n_particles_species;
  float   *x_temp;
  float   *y_temp;
  float   *z_temp;
  float   *h_smooth_temp;
  double   x_hat;
  double   y_hat;
  double   z_hat;
  double   theta;
  double   theta_roll;
  double   d_x_o;
  double   d_y_o;
  double   d_z_o;
  double   d_o;
  double   d_hat;
  double   particle_radius;
  double   x_tmp,y_tmp,z_tmp;
  int      flag_log;
  double   z_test;
  int      flag_use_Gadget;
  float    box_size_float=(float)box_size;
  float    half_box      =0.5*box_size_float;
  
  SID_log("Initializing projection-space...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Plane parallel projection?
  int flag_plane_parallel;
  if(check_mode_for_flag(camera_mode,CAMERA_PLANE_PARALLEL))
     flag_plane_parallel=TRUE;
  else
     flag_plane_parallel=FALSE;

  // Report the state of some flags
  if(flag_plane_parallel) SID_log("Plane-parallel  projection: ON", SID_LOG_COMMENT);
  else                    SID_log("Plane-parallel  projection: OFF",SID_LOG_COMMENT);
  if(flag_comoving)       SID_log("Comoving        projection: ON", SID_LOG_COMMENT);
  else                    SID_log("Comoving        projection: OFF",SID_LOG_COMMENT);
  if(flag_force_periodic) SID_log("Centre periodic projection: ON", SID_LOG_COMMENT);
  else                    SID_log("Centre periodic projection: OFF",SID_LOG_COMMENT);

  // Initialize the mapping quantities
  map_quantities_info mq;
  init_particle_map_quantities(&mq,render,render->camera->transfer_list,flag_comoving,expansion_factor);
  (*n_particles)       =mq.n_particles;
  (*flag_weigh)        =mq.flag_weigh;
  (*flag_line_integral)=mq.flag_line_integral;
  ptype_used           =mq.ptype_used;

  double x_o_in=x_o;
  double y_o_in=y_o;
  double z_o_in=z_o;
  double x_c_in=x_c;
  double y_c_in=y_c;
  double z_c_in=z_c;
  double FOV_x =FOV_x_in;
  double FOV_y =FOV_y_in;
  compute_perspective_transformation(x_o_in,
                                     y_o_in,
                                     z_o_in,
                                     x_c_in,
                                     y_c_in,
                                     z_c_in,
                                     unit_factor,
                                     unit_text,
                                     f_image_plane,
                                     stereo_offset,
                                     &FOV_x,
                                     &FOV_y,
                                     &d_o,
                                     &x_o,
                                     &y_o,
                                     &z_o,
                                     &x_c,
                                     &y_c,
                                     &z_c,
                                     &x_hat,
                                     &y_hat,
                                     &z_hat,
                                     &theta,
                                     &theta_roll);

  // The previous line sets d_o to the object distance.  Compute the distance to the image plane.
  double d_image_plane=d_o*f_image_plane;

  // Set mark arrays
  int    flag_mark_on=FALSE;
  char **mark        =(char **)SID_malloc(sizeof(char *)*N_GADGET_TYPE);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
      if(ADaPS_exist(render->plist_list[0]->data,"mark_%s",render->plist_list[0]->species[i_type])){
         mark[i_type]=(char *)ADaPS_fetch(render->plist_list[0]->data,"mark_%s",render->plist_list[0]->species[i_type]);
         flag_mark_on=TRUE;
      }
      else
         mark[i_type]=NULL;
  }

  // Determine how many particles are contributing to each 
  //    column of the image
  float   x_i;
  float   y_i;
  float   z_i;
  float   h_i;
  float   f_i;
  float   v_i;
  float   w_i;
  char    c_i;
  int     i_x;
  size_t *n_column;
  n_column=(size_t *)SID_calloc(sizeof(size_t)*nx);
  SID_log("Compute slab domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
    if(ptype_used[i_type] && ADaPS_exist(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type])){
      n_particles_species=((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type]))[0];
      for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){
         if(check_if_particle_marked(mark,i_type,i_particle,&c_i)){
            // Set the preoperties of the particle to be mapped (mode is TRUE because we need the angular size of the particle)
            set_particle_map_quantities(render,&mq,TRUE,k_particle,box_size_float,half_box,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

            // Transform particle to render-coordinates
            transform_particle(&x_i,
                               &y_i,
                               &z_i,
                               x_o,
                               y_o,
                               z_o,
                               x_hat,
                               y_hat,
                               z_hat,
                               d_o,
                               stereo_offset,
                               theta,
                               theta_roll,
                               box_size,
                               expansion_factor,
                               focus_shift_x,
                               focus_shift_y,
                               flag_comoving,
                               flag_force_periodic);

            if(z_i>d_near_field){
               double radius2_norm;
               double radius_kernel;
               double part_pos_x;
               double part_pos_y;
               int    kx_min;
               int    kx_max;
               int    ky_min;
               int    ky_max;
               int    n_ky;
               f_i=(float)compute_f_stretch(d_image_plane,z_i,flag_plane_parallel);
               set_pixel_space(h_i,
                               x_i,
                               y_i,
                               f_i,
                               xmin,
                               ymin,
                               FOV_x,
                               FOV_y,
                               pixel_size_x,
                               pixel_size_y,
                               radius_kernel_max,
                               &radius2_norm,
                               &radius_kernel,
                               &part_pos_x,
                               &part_pos_y,
                               &kx_min,
                               &kx_max,
                               &ky_min,
                               &ky_max);
               kx_min=MAX(kx_min,0);
               kx_max=MIN(kx_max,nx-1);
               ky_min=MAX(ky_min,0);
               ky_max=MIN(ky_max,ny-1);
               n_ky  =(ky_max-ky_min+1); 
               // It's not clear which of these two scalings should be used.  Scaling
               //   with no. of pixels yields imperfect results.
               if(n_ky>0){
                  for(i_x=kx_min;i_x<=kx_max;i_x++)
                     //n_column[i_x]+=n_ky; // Scale with no. of pixels
                     n_column[i_x]++;       // Scale with no. of particles
               }
            }
         }
      }
      j_particle+=n_particles_species;
    }
  }
  SID_Allreduce(SID_IN_PLACE,n_column,nx,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);

  // Create cumulative histogram
  for(i_x=1;i_x<nx;i_x++)
     n_column[i_x]+=n_column[i_x-1];

  // Decide on a range of columns to assign to each rank
  int    i_x_min_local;
  int    i_x_max_local;
  size_t norm     =n_column[nx-1];
  size_t n_used   =0;
  int    i_x_start=1;
  int    i_rank;
  size_t n_target;
  for(i_rank=0,i_x=0;i_rank<SID.n_proc;i_rank++){
     i_x_start=i_x;
     n_target =n_used+(size_t)((float)(norm-n_used)/(float)(SID.n_proc-i_rank));
     while(i_x<(nx-1) && n_column[i_x]<n_target) i_x++;
     if(i_rank==SID.My_rank){
        i_x_min_local=i_x_start;
        i_x_max_local=i_x;
     }
     n_used=n_column[i_x];
     i_x++;
  }
  if(SID.I_am_last_rank)
     if(i_x_max_local<nx)
        i_x_max_local=nx-1;
  SID_log("Done.",SID_LOG_CLOSE);

  // Create a couple arrays storing the domain decomposition for all
  int  i_x_min_send;
  int  i_x_max_send;
  int *i_x_min_rank;
  int *i_x_max_rank;
  i_x_min_rank=(int *)SID_malloc(sizeof(int)*SID.n_proc);
  i_x_max_rank=(int *)SID_malloc(sizeof(int)*SID.n_proc);
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
     i_x_min_send=i_x_min_local;
     i_x_max_send=i_x_max_local;
     SID_Bcast(&i_x_min_send,sizeof(int),i_rank,SID.COMM_WORLD);
     SID_Bcast(&i_x_max_send,sizeof(int),i_rank,SID.COMM_WORLD);
     i_x_min_rank[i_rank]=i_x_min_send;
     i_x_max_rank[i_rank]=i_x_max_send;
  }

  // Count the number of particles contributing to each rank
  size_t *n_rank;
  size_t *n_rank_local;
  size_t  n_particles_local;
  n_rank      =(size_t *)SID_calloc(sizeof(size_t)*SID.n_proc);
  n_rank_local=(size_t *)SID_calloc(sizeof(size_t)*SID.n_proc);
  SID_log("Count particles...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
    if(ptype_used[i_type] && ADaPS_exist(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type])){
      n_particles_species=((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type]))[0];
      for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){
         if(check_if_particle_marked(mark,i_type,i_particle,&c_i)){
            // Set the preoperties of the particle to be mapped
            set_particle_map_quantities(render,&mq,TRUE,k_particle,box_size_float,half_box,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

            // Transform particle to render-coordinates
            transform_particle(&x_i,
                               &y_i,
                               &z_i,
                               x_o,
                               y_o,
                               z_o,
                               x_hat,
                               y_hat,
                               z_hat,
                               d_o,
                               stereo_offset,
                               theta,
                               theta_roll,
                               box_size,
                               expansion_factor,
                               focus_shift_x,
                               focus_shift_y,
                               flag_comoving,
                               flag_force_periodic);

            if(z_i>d_near_field){
               double radius2_norm;
               double radius_kernel;
               double part_pos_x;
               double part_pos_y;
               int    kx_min;
               int    kx_max;
               int    ky_min;
               int    ky_max;
               int    n_ky;
               f_i=(float)compute_f_stretch(d_image_plane,z_i,flag_plane_parallel);
               set_pixel_space(h_i,
                               x_i,
                               y_i,
                               f_i,
                               xmin,
                               ymin,
                               FOV_x,
                               FOV_y,
                               pixel_size_x,
                               pixel_size_y,
                               radius_kernel_max,
                               &radius2_norm,
                               &radius_kernel,
                               &part_pos_x,
                               &part_pos_y,
                               &kx_min,
                               &kx_max,
                               &ky_min,
                               &ky_max);
               kx_min=MAX(kx_min,0);
               kx_max=MIN(kx_max,nx-1);
               ky_min=MAX(ky_min,0);
               ky_max=MIN(ky_max,ny-1);
               n_ky  =(ky_max-ky_min+1);
               if(n_ky>0){
                  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
                     if(kx_min<=i_x_max_rank[i_rank] && kx_max>=i_x_min_rank[i_rank])
                        n_rank_local[i_rank]++;
                  }
               }
            }
         }
      }
      j_particle+=n_particles_species;
    }
  }
  SID_Allreduce(n_rank_local,n_rank,SID.n_proc,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  n_particles_local=n_rank[SID.My_rank];
  (*n_particles)   =n_particles_local;
  SID_free(SID_FARG n_column);
  SID_log("Done.",SID_LOG_CLOSE);

  // Report decomposition results
  if(SID.n_proc>1){
     SID_log("Results of domain decomposition:",SID_LOG_OPEN);
     for(i_rank=0;i_rank<SID.n_proc;i_rank++)
        SID_log("Rank %03d image range: %5d->%5d n_particles=%zd",SID_LOG_COMMENT,i_rank,i_x_min_rank[i_rank],i_x_max_rank[i_rank],n_rank[i_rank]);
     SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
  }

  // Determine the needed size of the comm buffers
  size_t n_buffer;
  size_t n_buffer_max;
  calc_max(n_rank_local,&n_buffer,SID.n_proc,SID_SIZE_T,CALC_MODE_DEFAULT);
  SID_Allreduce(&n_buffer,&n_buffer_max,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);

  // Broadcast particles
  SID_log("Broadcast particles...(buffer=%zd particles)...",SID_LOG_OPEN|SID_LOG_TIMER,n_buffer_max);
  float  *x_buffer;
  float  *y_buffer;
  float  *z_buffer;
  float  *h_buffer;
  float  *f_buffer;
  float  *v_buffer;
  float  *w_buffer;
  char   *c_buffer;
  (*x)        =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*y)        =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*z)        =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*h_smooth) =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*f_stretch)=(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*value)    =(float *)SID_malloc(sizeof(float)*n_particles_local);
  (*weight)   =(float *)SID_malloc(sizeof(float)*n_particles_local);
  x_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  y_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  z_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  h_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  f_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  v_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  w_buffer    =(float *)SID_malloc(sizeof(float)*n_buffer);
  if(flag_mark_on){
     (*colour)=(char  *)SID_malloc(sizeof(char)*n_particles_local);
     c_buffer =(char  *)SID_malloc(sizeof(char)*n_buffer);
  }
  else{
     (*colour)=NULL;
     c_buffer =NULL;
  }
  size_t i_particle_rank;
  size_t j_particle_rank=0;
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
     int rank_to;
     int rank_from;
     set_exchange_ring_ranks(&rank_to,&rank_from,i_rank);
     i_particle_rank=0;
     for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
       if(ptype_used[i_type]){
         n_particles_species=((size_t *)ADaPS_fetch(render->plist_list[0]->data,"n_%s",render->plist_list[0]->species[i_type]))[0];
         for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){
            if(check_if_particle_marked(mark,i_type,i_particle,&c_i)){
               // Set the preoperties of the particle to be mapped
               set_particle_map_quantities(render,&mq,TRUE,k_particle,box_size_float,half_box,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

               // Transform particle to render-coordinates
               transform_particle(&x_i,
                                  &y_i,
                                  &z_i,
                                  x_o,
                                  y_o,
                                  z_o,
                                  x_hat,
                                  y_hat,
                                  z_hat,
                                  d_o,
                                  stereo_offset,
                                  theta,
                                  theta_roll,
                                  box_size,
                                  expansion_factor,
                                  focus_shift_x,
                                  focus_shift_y,
                                  flag_comoving,
                                  flag_force_periodic);

               if(z_i>d_near_field){
                  double radius2_norm;
                  double radius_kernel;
                  double part_pos_x;
                  double part_pos_y;
                  int    kx_min;
                  int    kx_max;
                  int    ky_min;
                  int    ky_max;
                  int    n_ky;
                  f_i=(float)compute_f_stretch(d_image_plane,z_i,flag_plane_parallel);
                  set_pixel_space(h_i,
                                  x_i,
                                  y_i,
                                  f_i,
                                  xmin,
                                  ymin,
                                  FOV_x,
                                  FOV_y,
                                  pixel_size_x,
                                  pixel_size_y,
                                  radius_kernel_max,
                                  &radius2_norm,
                                  &radius_kernel,
                                  &part_pos_x,
                                  &part_pos_y,
                                  &kx_min,
                                  &kx_max,
                                  &ky_min,
                                  &ky_max);
                  kx_min=MAX(kx_min,0);
                  kx_max=MIN(kx_max,nx-1);
                  ky_min=MAX(ky_min,0);
                  ky_max=MIN(ky_max,ny-1);
                  n_ky  =(ky_max-ky_min+1);
                  if(n_ky>0){
                     if(kx_min<=i_x_max_rank[rank_to] && kx_max>=i_x_min_rank[rank_to]){
                        x_buffer[i_particle_rank]=x_i;
                        y_buffer[i_particle_rank]=y_i;
                        z_buffer[i_particle_rank]=z_i;
                        h_buffer[i_particle_rank]=h_i;
                        f_buffer[i_particle_rank]=f_i;
                        w_buffer[i_particle_rank]=w_i;
                        v_buffer[i_particle_rank]=v_i;
                        if(c_buffer!=NULL)
                           c_buffer[i_particle_rank]=(char)c_i;
                        i_particle_rank++;
                     }
                  }
               }
            }
         }
         j_particle+=n_particles_species;
       }
     } // loop over particles

     // Perform exchanges
     size_t n_exchange;
     exchange_ring_buffer(x_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*x)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(y_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*y)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(z_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*z)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(h_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*h_smooth)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(f_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*f_stretch)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(w_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*weight)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     exchange_ring_buffer(v_buffer,
                          sizeof(float),
                          n_rank_local[rank_to],
                          &((*value)[j_particle_rank]),
                          &n_exchange,
                          i_rank);
     if(c_buffer!=NULL)
        exchange_ring_buffer(c_buffer,
                             sizeof(char),
                             n_rank_local[rank_to],
                             &((*colour)[j_particle_rank]),
                             &n_exchange,
                             i_rank);
     j_particle_rank+=n_exchange;
  }

  // Sanity checks
  if(j_particle_rank!=n_particles_local)
     SID_trap_error("The wrong number of particles were received (ie. %zd!=%zd) on rank %d.",ERROR_LOGIC,
                    j_particle_rank,n_particles_local,SID.My_rank);

  SID_log("Done.",SID_LOG_CLOSE);

  // Sort the local particles by position
  merge_sort((*z),(size_t)(*n_particles),z_index,SID_FLOAT,SORT_COMPUTE_INDEX,FALSE);

  // Clean-up
  SID_free(SID_FARG n_rank);
  SID_free(SID_FARG n_rank_local);
  SID_free(SID_FARG i_x_min_rank);
  SID_free(SID_FARG i_x_max_rank);
  SID_free(SID_FARG x_buffer);
  SID_free(SID_FARG y_buffer);
  SID_free(SID_FARG z_buffer);
  SID_free(SID_FARG h_buffer);
  SID_free(SID_FARG f_buffer);
  SID_free(SID_FARG v_buffer);
  SID_free(SID_FARG w_buffer);
  SID_free(SID_FARG c_buffer);
  free_particle_map_quantities(&mq);
  SID_free(SID_FARG mark);

  (*i_x_min_local_return)=i_x_min_local;
  (*i_x_max_local_return)=i_x_max_local;

  SID_log("Done.",SID_LOG_CLOSE);
}

void fetch_image_array(image_info *image,double **values);
void fetch_image_array(image_info *image,double **values){
   if(image!=NULL)
      (*values)=image->values;
   else
      (*values)=NULL;
}

void render_frame(render_info  *render){
  size_t     i_particle;
  size_t     j_particle;
  size_t     k_particle;
  int        i_pixel;              
  size_t     i_kernel;
  size_t     n_pixels;
  size_t     n_unmasked;
  double     xmin,ymin;
  double     pixel_size;
  double     pixel_size_x;
  double     pixel_size_y;
  double     pixel_area;
  float     *weight =NULL;
  float     *value  =NULL;
  char      *colour =NULL;
  size_t    *z_index=NULL;
  int        kx_min;
  int        kx_max;
  int        ky_min;
  int        ky_max;
  int        kz_min;
  int        kz_max;
  char      *mask       =NULL;
  double    *temp_image =NULL;
  double    *RGB_image  =NULL;
  double    *Y_image    =NULL;
  double    *z_image    =NULL;
  double    *RY_image   =NULL;
  double    *GY_image   =NULL;
  double    *BY_image   =NULL;
  double     part_pos_x;
  double     part_pos_y;
  double     z_i;
  double     part_h_z;
  double     part_h_xy;
  double     pixel_pos_x;
  double     pixel_pos_y;
  double     prob;
  int        kx;
  int        ky;
  int        pos;
  size_t     n_particles;
  float      *x;
  float      *y;
  float      *z;
  float      *h_smooth;
  float      *f_stretch;
  size_t     n_x;
  size_t     n_y;
  size_t     n_z;
  double     kernel;
  int        flag_weigh;
  int        flag_line_integral;
  int        n_images;

  double     deldr2;
  double     deldr2i;
  double     radius2;
  double     radius4;
  double     radius;
  double     ikernel;
  double     d2;
  double     d;
  double     radius_kernel;
  double     radius2_norm;
  double     f_table;
  int        i_table;
  double    *kernel_radius;
  double    *kernel_table;
  double     kernel_table_avg;
  double     radius_kernel_max;
  double     radius_kernel_max2;

  double     d_o,d_x_o,d_y_o,d_z_o;

  double       x_o;
  double       y_o;
  double       z_o;
  double       x_c;
  double       y_c;
  double       z_c;
  double       FOV;
  double       FOV_x_object_plane;
  double       FOV_y_object_plane;
  double       FOV_x_image_plane;
  double       FOV_y_image_plane;
  double       box_size;
  int          nx;
  int          ny;
  int          v_mode;
  int          w_mode;
  char        *parameter;
  int          flag_comoving;
  int          flag_fade;
  double       alpha_fade;
  int          flag_force_periodic;
  int          flag_add_absorption;
  double       expansion_factor;
  double       focus_shift_x;
  double       focus_shift_y;
  ADaPS       *transfer;
  double       h_Hubble;
  double       f_image_plane;
  double       d_image_plane;
  double       d_near_field;
  double       d_taper_field;
  double       stereo_offset;
  int          i_x,i_y;
  int         *mask_buffer;

  int          i_image;
  int          camera_mode;

  double       f_absorption;
  f_absorption=render->f_absorption;
  if(f_absorption<0.)
     f_absorption=0.;

  x_o          =render->camera->perspective->p_o[0];
  y_o          =render->camera->perspective->p_o[1];
  z_o          =render->camera->perspective->p_o[2];
  x_c          =render->camera->perspective->p_c[0];
  y_c          =render->camera->perspective->p_c[1];
  z_c          =render->camera->perspective->p_c[2];
  d_o          =render->camera->perspective->d_o;
  FOV          =render->camera->perspective->FOV;
  focus_shift_x=render->camera->perspective->focus_shift_x;
  focus_shift_y=render->camera->perspective->focus_shift_y;
  h_Hubble     =render->h_Hubble;

  double unit_factor;
  char   unit_text[32];
  int    kernel_flag;
  int    flag_velocity_space=render->camera->flag_velocity_space;
  if(flag_velocity_space){
     SID_log("Rendering in velocity-space.",SID_LOG_COMMENT);
     kernel_flag=SPH_KERNEL_GAUSSIAN;
     unit_factor=1e3;
     sprintf(unit_text,"km/s");
  }
  else{
     SID_log("Rendering in real-space.",SID_LOG_COMMENT);
     kernel_flag=SPH_KERNEL_GADGET;
     unit_factor=M_PER_MPC/h_Hubble;
     sprintf(unit_text,"Mpc/h");
  }

  // Put all the render distances into the Gadget units of the
  //    relevant space (r or v-space)
  x_o          *=unit_factor;
  y_o          *=unit_factor;
  z_o          *=unit_factor;
  x_c          *=unit_factor;
  y_c          *=unit_factor;
  z_c          *=unit_factor;
  d_o          *=unit_factor;
  FOV          *=unit_factor;
  focus_shift_x*=unit_factor;
  focus_shift_y*=unit_factor;

  nx                 =render->camera->width;
  ny                 =render->camera->height;
  v_mode             =render->v_mode;
  w_mode             =render->w_mode;
  camera_mode        =render->camera->camera_mode;
  flag_comoving      =render->flag_comoving;
  flag_fade          =render->flag_fade;
  alpha_fade         =render->alpha_fade;
  flag_force_periodic=render->flag_force_periodic;
  flag_add_absorption=render->flag_add_absorption;
  expansion_factor   =render->camera->perspective->time;
  f_image_plane      =render->camera->f_image_plane;
  d_near_field       =render->camera->f_near_field*d_o;
  d_taper_field      =render->camera->f_taper_field*d_o;
  d_image_plane      =d_o*f_image_plane;
  box_size           =((double *)ADaPS_fetch(render->plist_list[0]->data,"box_size"))[0];

  // Sanity check the taper distances
  double taper_width=d_taper_field-d_near_field;
  if(taper_width<0.)
     SID_trap_error("The near-field distance (%le [%s]) must be less than the taper distance (%le [%s]) if non-zero.",ERROR_LOGIC,
                    d_near_field*unit_factor, unit_text,
                    d_taper_field*unit_factor,unit_text);

  // Set FOV
  if(d_near_field>0.)
     SID_log("Near  field  = %le [%s]",SID_LOG_COMMENT,d_near_field *unit_factor,unit_text);
  if(d_taper_field>0.)
     SID_log("Taper field  = %le [%s]",SID_LOG_COMMENT,d_taper_field*unit_factor,unit_text);
  if(d_near_field>0. || d_taper_field>0.)
     SID_log("Image plane  = %le [%s]",SID_LOG_COMMENT,d_image_plane*unit_factor,unit_text);
  SID_log("f_absorption = %le",SID_LOG_COMMENT,f_absorption);
  if(nx>=ny){
    FOV_y_object_plane=FOV;
    FOV_x_object_plane=FOV_y_object_plane*(double)nx/(double)ny;    
  }
  else{
    FOV_x_object_plane=FOV;
    FOV_y_object_plane=FOV_x_object_plane*(double)nx/(double)ny;        
  }
  FOV_x_image_plane=FOV_x_object_plane*f_image_plane;
  FOV_y_image_plane=FOV_y_object_plane*f_image_plane;

  // Loop over the left/right stereo pair (if necessary)
  if(check_mode_for_flag(camera_mode,CAMERA_STEREO))
    i_image =0;
  else
    i_image =1;
  for(;i_image<2;i_image++){

    switch(i_image){
      // Left image
      case 0:
        fetch_image_array(render->camera->image_RGB_left, &RGB_image);
        fetch_image_array(render->camera->image_RGBY_left,&temp_image); // use as workspace
        fetch_image_array(render->camera->image_Y_left,   &Y_image);
        fetch_image_array(render->camera->image_Z_left,   &z_image);
        fetch_image_array(render->camera->image_RY_left,  &RY_image);
        fetch_image_array(render->camera->image_GY_left,  &GY_image);
        fetch_image_array(render->camera->image_BY_left,  &BY_image);
        stereo_offset=-d_image_plane/render->camera->stereo_ratio;
        break;
      case 1:
        // Right image
        if(check_mode_for_flag(camera_mode,CAMERA_STEREO)){
           fetch_image_array(render->camera->image_RGB_right, &RGB_image);
           fetch_image_array(render->camera->image_RGBY_right,&temp_image); // use as workspace
           fetch_image_array(render->camera->image_Y_right,   &Y_image);
           fetch_image_array(render->camera->image_Z_right,   &z_image);
           fetch_image_array(render->camera->image_RY_right,  &RY_image);
           fetch_image_array(render->camera->image_GY_right,  &GY_image);
           fetch_image_array(render->camera->image_BY_right,  &BY_image);
           stereo_offset=d_image_plane/render->camera->stereo_ratio;
        }
        // Mono image (stereo turned off)
        else{
           fetch_image_array(render->camera->image_RGB, &RGB_image);
           fetch_image_array(render->camera->image_RGBY,&temp_image); // use as workspace
           fetch_image_array(render->camera->image_Y,   &Y_image);
           fetch_image_array(render->camera->image_Z,   &z_image);
           fetch_image_array(render->camera->image_RY,  &RY_image);
           fetch_image_array(render->camera->image_GY,  &GY_image);
           fetch_image_array(render->camera->image_BY,  &BY_image);
           stereo_offset=0;           
        }
        break;
    }

    SID_log("Projecting to a %dx%d pixel array...",SID_LOG_OPEN|SID_LOG_TIMER,nx,ny);

    // Set physical image-plane domain
    xmin  = -FOV_x_image_plane/2.; // Things will be centred on (x_o,y_o,z_o) later
    ymin  = -FOV_y_image_plane/2.; // Things will be centred on (x_o,y_o,z_o) later

    xmin-=stereo_offset;

    // Compute image scales
    n_pixels    =nx*ny;
    pixel_size_x=FOV_x_image_plane/(double)nx;
    pixel_size_y=FOV_y_image_plane/(double)ny;
    pixel_size  =0.5*(pixel_size_x+pixel_size_y);
    pixel_area  =pixel_size_x*pixel_size_y;
    if(fabs((pixel_size_x-pixel_size_y)/pixel_size_x)>1e-4)
      SID_log_warning("pixels are not square by %7.3f%%",0,fabs((pixel_size_x-pixel_size_y)/pixel_size_x)*1e2);

    // Generate the smoothing kernal
    set_sph_kernel(&(render->kernel_radius),
                   &(render->kernel_table_3d),
                   &(render->kernel_table),
                   &(render->kernel_table_avg),
                   kernel_flag|SPH_KERNEL_2D);
    kernel_radius     =render->kernel_radius;
    kernel_table      =render->kernel_table;
    kernel_table_avg  =render->kernel_table_avg;
    radius_kernel_max =kernel_radius[N_KERNEL_TABLE];
    radius_kernel_max2=radius_kernel_max*radius_kernel_max;

    double x_o_in=x_o;
    double y_o_in=y_o;
    double z_o_in=z_o;
    double x_c_in=x_c;
    double y_c_in=y_c;
    double z_c_in=z_c;
    double FOV_x =FOV_x_object_plane;
    double FOV_y =FOV_y_object_plane;
    double x_o_out;
    double y_o_out;
    double z_o_out;
    double x_c_out;
    double y_c_out;
    double z_c_out;
    double x_hat;
    double y_hat;
    double z_hat;
    double theta;
    double theta_roll;
    double d_o;
    compute_perspective_transformation(x_o_in,
                                       y_o_in,
                                       z_o_in,
                                       x_c_in,
                                       y_c_in,
                                       z_c_in,
                                       unit_factor,
                                       unit_text,
                                       f_image_plane,
                                       stereo_offset,
                                       &FOV_x,
                                       &FOV_y,
                                       &d_o, // Includes f_image_plane correction
                                       &x_o_out,
                                       &y_o_out,
                                       &z_o_out,
                                       &x_c_out,
                                       &y_c_out,
                                       &z_c_out,
                                       &x_hat,
                                       &y_hat,
                                       &z_hat,
                                       &theta,
                                       &theta_roll);

    // Initialize make_map
    int i_x_min_local;
    int i_x_max_local;
    if(flag_add_absorption)
       init_make_map_abs(render,
                         x_o,y_o,z_o,
                         x_c,y_c,z_c,
                         unit_factor,unit_text,
                         f_image_plane,
                         box_size,FOV_x_object_plane,FOV_y_object_plane,
                         xmin,ymin,
                         pixel_size_x,pixel_size_y,
                         radius_kernel_max,
                         nx,ny,
                         expansion_factor,
                         focus_shift_x,
                         focus_shift_y,
                         d_near_field,
                         stereo_offset,
                         flag_comoving,
                         flag_force_periodic,
                         camera_mode,
                         &flag_weigh,
                         &flag_line_integral,
                         &x,&y,&z,
                         &h_smooth,
                         &f_stretch,
                         &value,
                         &weight,
                         &colour,
                         &z_index,
                         &i_x_min_local,
                         &i_x_max_local,
                         &n_particles);
    else
       init_make_map_noabs(render,
                           x_o,y_o,z_o,
                           x_c,y_c,z_c,
                           unit_factor,unit_text,
                           f_image_plane,
                           box_size,FOV_x_object_plane,FOV_y_object_plane,
                           xmin,ymin,
                           pixel_size_x,pixel_size_y,
                           radius_kernel_max,
                           nx,ny,
                           expansion_factor,
                           focus_shift_x,
                           focus_shift_y,
                           d_near_field,
                           stereo_offset,
                           flag_comoving,
                           flag_force_periodic,
                           camera_mode,
                           &flag_weigh,
                           &flag_line_integral,
                           &x,&y,&z,
                           &h_smooth,
                           &f_stretch,
                           &value,
                           &weight,
                           &colour,
                           &z_index,
                           &i_x_min_local,
                           &i_x_max_local,
                           &n_particles);

    // Initialize image arrays
    mask=(char *)SID_malloc(sizeof(char)*n_pixels);
    for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
      if(RGB_image!=NULL)
         RGB_image[i_pixel]=0.;
      if(temp_image!=NULL)
         temp_image[i_pixel]=0.;
      if(Y_image!=NULL)
         Y_image[i_pixel]=0.;
      if(z_image!=NULL)
        z_image[i_pixel]=0.;
      if(RY_image!=NULL)
        RY_image[i_pixel]=0.;
      if(GY_image!=NULL)
        GY_image[i_pixel]=0.;
      if(BY_image!=NULL)
        BY_image[i_pixel]=0.;
      mask[i_pixel]=FALSE;
    }

    // Perform projection
    size_t        n_particles_used_local=0;
    size_t        n_particles_used      =0;
    pcounter_info pcounter;
    SID_init_pcounter(&pcounter,n_particles,10);
    SID_log("Performing projection...",SID_LOG_OPEN|SID_LOG_TIMER);
    n_particles_used_local=0;
    for(size_t ii_particle=0;ii_particle<n_particles;ii_particle++){
      i_particle=z_index[n_particles-1-ii_particle];
      z_i       =(double)z[i_particle];

      part_h_z=(double)h_smooth[i_particle];

      if(z_i>d_near_field){
        n_particles_used_local++;

        // Set pixel space ranges and positions
        part_h_xy    =part_h_z*f_stretch[i_particle];
        radius2_norm =1./(part_h_xy*part_h_xy);
        part_pos_x   =(double)(x[i_particle]*f_stretch[i_particle]);
        part_pos_y   =(double)(y[i_particle]*f_stretch[i_particle]);
        radius_kernel=part_h_xy;
        kx_min=(int)((part_pos_x-radius_kernel-xmin)/pixel_size_x);
        kx_max=(int)((part_pos_x+radius_kernel-xmin)/pixel_size_x+ONE_HALF);
        ky_min=(int)((part_pos_y-radius_kernel-ymin)/pixel_size_y);
        ky_max=(int)((part_pos_y+radius_kernel-ymin)/pixel_size_y+ONE_HALF);

        // Compute any potential fading
        // alpha_fade=2 for normal inverse-square fading past the image plane
        double f_fade;
        if(flag_fade && z_i>d_o)
           f_fade=pow(d_image_plane/z_i,-alpha_fade);
        else
           f_fade=1;

        // Compute any potential tapering
        double f_taper;
        if(taper_width>0. && z_i<d_taper_field)
           f_taper=(z_i-d_near_field)/taper_width;
        else
           f_taper=1;

        // Combine dimming factors into one
        double f_dim;
        f_dim=f_taper*f_fade;

        // Set the particle values and weights
        double v_i =(double)value[i_particle];
        double w_i =(double)weight[i_particle];
        double R_i =1.;
        double G_i =1.;
        double B_i =1.;
        if(colour!=NULL){
           R_i =RGB_lookup(render,colour[i_particle],0);
           G_i =RGB_lookup(render,colour[i_particle],1);
           B_i =RGB_lookup(render,colour[i_particle],2);
        }
        double vw_i=v_i*w_i;
        double zw_i=z_i*w_i;

        // Loop over the kernal
//printf("%le %d\n",
//printf("%f %f %f -- %le %le -- %le %le -- %le %le %le -- %4d %4d %4d %4d -- %le %le %le\n",
//       part_pos_x,part_pos_y,f_stretch[i_particle],
//       radius_kernel,f_table,
//       pixel_size_x,pixel_size_y,
//       xmin,ymin,FOV,
//       kx_min,kx_max-kx_min,ky_min,ky_max-ky_min,
//       w_i,v_i,kernel);
//if(n_particles_used_local>=10) SID_exit(ERROR_NONE);
        for(kx=kx_min,pixel_pos_x=xmin+(kx_min+0.5)*pixel_size_x;kx<=kx_max;kx++,pixel_pos_x+=pixel_size_x){
          if(kx>=0 && kx<nx){
            for(ky=ky_min,pixel_pos_y=ymin+(ky_min+0.5)*pixel_size_y;ky<=ky_max;ky++,pixel_pos_y+=pixel_size_y){
              if(ky>=0 && ky<ny){
                radius2=
                  (pixel_pos_x-part_pos_x)*(pixel_pos_x-part_pos_x)+
                  (pixel_pos_y-part_pos_y)*(pixel_pos_y-part_pos_y);
                radius2*=radius2_norm;
                // Construct images here
                if(radius2<radius_kernel_max2){
                  pos    =ky+kx*ny;
                  f_table=sqrt(radius2)/radius_kernel_max;
                  i_table=(int)(f_table*(double)N_KERNEL_TABLE);
                  kernel =kernel_table[i_table]+
                    (kernel_table[i_table+1]-kernel_table[i_table])*
                    (f_table-kernel_radius[i_table])*(double)N_KERNEL_TABLE;
                  double w_k=w_i*kernel;
                  // Perform addition
                  if(Y_image!=NULL)    Y_image[pos]   +=(f_dim    -f_absorption*   Y_image[pos])*w_k;
                  if(temp_image!=NULL) temp_image[pos]+=(f_dim*v_i-f_absorption*temp_image[pos])*w_k;
                  if(z_image!=NULL)    z_image[pos]   +=(f_dim*z_i-f_absorption*   z_image[pos])*w_k;
                  if(RY_image!=NULL)   RY_image[pos]  +=(f_dim*R_i-f_absorption*  RY_image[pos])*w_k;
                  if(GY_image!=NULL)   GY_image[pos]  +=(f_dim*G_i-f_absorption*  GY_image[pos])*w_k;
                  if(BY_image!=NULL)   BY_image[pos]  +=(f_dim*B_i-f_absorption*  BY_image[pos])*w_k;
                  mask[pos]=TRUE;
                }
              }
            } // loop over kernel
          } 
        } // loop over kernel
      } // near-field selection
      SID_check_pcounter(&pcounter,ii_particle);
    } // loop over particles
    SID_Barrier(SID.COMM_WORLD);
    SID_Allreduce(&n_particles_used_local,&n_particles_used,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
    SID_log("n_particles_used=%zd",SID_LOG_COMMENT,n_particles_used);
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Image normalization, etc...",SID_LOG_OPEN|SID_LOG_TIMER);

    // Add results from all ranks if this is being run in parallel
    //   First, clear the parts of the image not in this rank's domain ...
    for(kx=0;kx<nx;kx++){
       if(!(kx>=i_x_min_local && kx<=i_x_max_local)){
          for(ky=0;ky<ny;ky++){
             pos=ky+kx*ny;
             if(temp_image!=NULL)
                temp_image[pos]=0.;
             if(Y_image!=NULL)
                Y_image[pos]=0.;
             if(z_image!=NULL)
                z_image[pos]=0.;
             if(RY_image!=NULL){
                RY_image[pos]=0.;
                GY_image[pos]=0.;
                BY_image[pos]=0.;
             }
          }
       }
    }

    // Join mask results.  The MPI 1.1 standard does not specify MPI_SUM 
    //   on MPI_CHAR types so we need to do this awkward thing for the mask array ...
#if USE_MPI
    mask_buffer=(int *)SID_malloc(sizeof(int)*nx);
    for(i_y=0;i_y<ny;i_y++){
      for(i_x=0,i_pixel=i_y*nx;i_x<nx;i_x++,i_pixel++)
        mask_buffer[i_x]=(int)mask[i_pixel];
      SID_Allreduce(SID_IN_PLACE,mask_buffer,nx,SID_INT,SID_MAX,SID.COMM_WORLD);
      for(i_x=0,i_pixel=i_y*nx;i_x<nx;i_x++,i_pixel++){
         if(mask_buffer[i_x])
           mask[i_pixel]=TRUE;
      }
    }
    SID_free(SID_FARG mask_buffer);
#endif
    SID_Allreduce(SID_IN_PLACE,temp_image,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
    SID_Allreduce(SID_IN_PLACE,Y_image,   n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
    if(z_image!=NULL)
      SID_Allreduce(SID_IN_PLACE,z_image,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
    if(RY_image!=NULL){
      SID_Allreduce(SID_IN_PLACE,RY_image,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,GY_image,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,BY_image,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
    }

    // Create final normalized images and clear the temp_image which has been used as a buffer
    if(RGB_image!=NULL){
       for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
             RGB_image[i_pixel]=temp_image[i_pixel]/Y_image[i_pixel];
          temp_image[i_pixel]=0.;
       }
    }
    if(RY_image!=NULL){
       for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
             RY_image[i_pixel]=RY_image[i_pixel]/Y_image[i_pixel];
       }
    }
    if(GY_image!=NULL){
       for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
             GY_image[i_pixel]=GY_image[i_pixel]/Y_image[i_pixel];
       }
    }
    if(BY_image!=NULL){
       for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
             BY_image[i_pixel]=BY_image[i_pixel]/Y_image[i_pixel];
       }
    }

    // Normalize z-frame (if needed)
    if(z_image!=NULL){
       for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
         if(mask[i_pixel])
           z_image[i_pixel]=z_image[i_pixel]/Y_image[i_pixel];
         else
           z_image[i_pixel]=0.;
       }
    }

    // Take log_10 (if needed)
    if(check_mode_for_flag(v_mode,MAKE_MAP_LOG)){
      SID_log("Taking log of RGB image...",SID_LOG_OPEN);
      for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
        if(mask[i_pixel])
          RGB_image[i_pixel]=take_log10(RGB_image[i_pixel]);
        else
          RGB_image[i_pixel]=LOG_ZERO;
      }
      SID_log("Done.",SID_LOG_CLOSE);
    }
    if(check_mode_for_flag(w_mode,MAKE_MAP_LOG)){
      SID_log("Taking log of Y image...",SID_LOG_OPEN);
      for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
        if(mask[i_pixel])
          Y_image[i_pixel]=take_log10(Y_image[i_pixel]);
        else
          Y_image[i_pixel]=LOG_ZERO;
      }
      SID_log("Done.",SID_LOG_CLOSE);
    }

    // Compute some image statistics
    double min_RGB_image;
    double max_RGB_image;
    double min_Y_image;
    double max_Y_image;
    double min_z_image;
    double max_z_image;
    double min_RY_image;
    double max_RY_image;
    double min_GY_image;
    double max_GY_image;
    double min_BY_image;
    double max_BY_image;
    for(i_pixel=0,n_unmasked=0;i_pixel<n_pixels;i_pixel++) 
      if(mask[i_pixel]) 
        n_unmasked++;
    if(RGB_image!=NULL){
      calc_min(RGB_image,&min_RGB_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(RGB_image,&max_RGB_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    if(Y_image!=NULL){
      calc_min(Y_image,&min_Y_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(Y_image,&max_Y_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    if(z_image!=NULL){
      calc_min(z_image,&min_z_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(z_image,&max_z_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    if(RY_image!=NULL){
      calc_min(RY_image,&min_RY_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(RY_image,&max_RY_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    if(GY_image!=NULL){
      calc_min(GY_image,&min_GY_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(GY_image,&max_GY_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    if(BY_image!=NULL){
      calc_min(BY_image,&min_BY_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(BY_image,&max_BY_image,n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  
    // Report image statistics
    if(n_unmasked>0){
       SID_log("Image statistics:",SID_LOG_OPEN);
       SID_log("coverage=%3d%%",SID_LOG_COMMENT,(int)(100.*(double)n_unmasked/(double)n_pixels));
       if(RGB_image!=NULL){
          SID_log("RGB min =%le",SID_LOG_COMMENT,min_RGB_image);
          SID_log("RGB max =%le",SID_LOG_COMMENT,max_RGB_image);
       }
       if(Y_image!=NULL){
          SID_log("Y   min =%le",SID_LOG_COMMENT,min_Y_image);
          SID_log("Y   max =%le",SID_LOG_COMMENT,max_Y_image);
       }
       if(RY_image!=NULL){
          SID_log("RY  min =%le",SID_LOG_COMMENT,min_RY_image);
          SID_log("RY  max =%le",SID_LOG_COMMENT,max_RY_image);
       }
       if(GY_image!=NULL){
          SID_log("GY  min =%le",SID_LOG_COMMENT,min_GY_image);
          SID_log("GY  max =%le",SID_LOG_COMMENT,max_GY_image);
       }
       if(BY_image!=NULL){
          SID_log("BY  min =%le",SID_LOG_COMMENT,min_BY_image);
          SID_log("BY  max =%le",SID_LOG_COMMENT,max_BY_image);
       }
       if(z_image!=NULL){
          SID_log("Z   min =%le [%s]",SID_LOG_COMMENT,min_z_image/unit_factor,unit_text);
          SID_log("Z   max =%le [%s]",SID_LOG_COMMENT,max_z_image/unit_factor,unit_text);
       }
       SID_log("",SID_LOG_SILENT_CLOSE);
    }
    else
      SID_log("IMAGE IS EMPTY.",SID_LOG_COMMENT);

    // Clean-up
    SID_free(SID_FARG x);
    SID_free(SID_FARG y);
    SID_free(SID_FARG z);
    SID_free(SID_FARG h_smooth);
    SID_free(SID_FARG f_stretch);
    SID_free(SID_FARG value);
    SID_free(SID_FARG weight);
    SID_free(SID_FARG colour);
    SID_free(SID_FARG z_index);
    SID_free(SID_FARG mask);
  
    SID_log("Done.",SID_LOG_CLOSE);
  }
}

