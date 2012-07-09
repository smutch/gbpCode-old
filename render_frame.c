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
                    double  radius_kernel_norm,
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
                    double  radius_kernel_norm,
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
   (*radius_kernel)=radius_kernel_norm*part_h_xy;
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
   
   // Shift image plane
   (*z_i)+=d_o;

}

#define MAKE_MAP_NO_WEIGHTING 0
#define MAKE_MAP_MODE_RHO     2
#define MAKE_MAP_MODE_SIGMA   4

typedef struct map_quantities_info map_quantities_info;
struct map_quantities_info{
   int           flag_weigh;
   int           flag_line_integral;
   int          *ptype_used;
   size_t        n_particles;
   GBPREAL      *x;
   GBPREAL      *y;
   GBPREAL      *z;
   float        *h_smooth;
   float        *rho;
   float        *sigma;
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
   SID_free(SID_FARG mq->ptype_used);
}

void init_particle_map_quantities(map_quantities_info *mq,plist_info *plist,ADaPS *transfer_list,char *parameter,int flag_comoving,double expansion_factor);
void init_particle_map_quantities(map_quantities_info *mq,plist_info *plist,ADaPS *transfer_list,char *parameter,int flag_comoving,double expansion_factor){

  SID_log("Initializing quantities...",SID_LOG_OPEN);

  // Defaults
  mq->flag_weigh             =FALSE;
  mq->flag_line_integral     =FALSE;
  mq->n_particles            =0;
  mq->h_smooth               =NULL;
  mq->x                      =NULL;
  mq->y                      =NULL;
  mq->z                      =NULL;
  mq->mass_array             =0.;
  mq->rho                    =NULL;
  mq->transfer_rho           =NULL;
  mq->flag_transfer_rho_log  =FALSE;
  mq->sigma                  =NULL;
  mq->transfer_sigma         =NULL;
  mq->flag_transfer_sigma_log=FALSE;

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

  // Render a dark matter surface density map
  if(!strcmp(parameter,"Sigma_M_dark")){
    mq->ptype_used[GADGET_TYPE_DARK]=TRUE;
    mq->flag_weigh                  =FALSE;
    mq->flag_line_integral          =FALSE;
    mq->n_particles                 =((size_t *)ADaPS_fetch(plist->data,"n_dark"))[0];
    if(ADaPS_exist(plist->data,"rho_dark"))
      mq->rho=(float *)ADaPS_fetch(plist->data,"rho_dark");
    else
      SID_trap_error("No denisities available to compute %s in make_map.",ERROR_LOGIC,parameter);
  }

  // Render a dark matter column depth map
  else if(!strcmp(parameter,"tau_dark")){
    mq->ptype_used[GADGET_TYPE_DARK]=TRUE;
    mq->flag_weigh                  =FALSE;
    mq->flag_line_integral          =TRUE;
    mq->n_particles                 =((size_t *)ADaPS_fetch(plist->data,"n_dark"))[0];
    mq->v_mode=MAKE_MAP_MODE_RHO;
    mq->w_mode=0;
    if(ADaPS_exist(plist->data,"rho_dark")){
      mq->rho=(float *)ADaPS_fetch(plist->data,"rho_dark");
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
    else if(ADaPS_exist(plist->data,"mass_array_dark")){
      mq->rho       =NULL;
      mq->mass_array=((double *)ADaPS_fetch(plist->data,"mass_array_dark"))[0];
    }
    else
      SID_trap_error("No densities or masses available to compute %s in make_map.",ERROR_LOGIC,parameter);
  }

  // Render a dark matter velocity dispersion map
  else if(!strcmp(parameter,"sigma_v_dark")){
    mq->ptype_used[GADGET_TYPE_DARK]=TRUE;
    mq->flag_weigh                  =TRUE;
    mq->flag_line_integral          =TRUE;
    mq->n_particles                 =((size_t *)ADaPS_fetch(plist->data,"n_dark"))[0];
    mq->v_mode=MAKE_MAP_MODE_SIGMA;
    mq->w_mode=MAKE_MAP_MODE_RHO;
    if(ADaPS_exist(plist->data,"rho_dark")){
      mq->rho  =(float *)ADaPS_fetch(plist->data,"rho_dark");
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
    else if(ADaPS_exist(plist->data,"mass_array_dark")){
      mq->rho       =NULL;
      mq->mass_array=((double *)ADaPS_fetch(plist->data,"mass_array_dark"))[0];
    }
    else
      SID_trap_error("No densities or masses available to compute %s in make_map.",ERROR_LOGIC,parameter);

    // Use sigma_v for values
    if(ADaPS_exist(plist->data,"sigma_v_dark")){
      mq->sigma=(float *)ADaPS_fetch(plist->data,"sigma_v_dark");
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
      SID_trap_error("No sigma_v's available to render %s in make_map.",ERROR_LOGIC,parameter);
  }
  else
    SID_trap_error("Unknown parameter {%s}.",ERROR_LOGIC,parameter);

  // Initialize smoothings
  int n_type_used=0;
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
     if(mq->ptype_used[i_type]){
        // This code can use one-and-only-one particle type at a time at the moment
        n_type_used++;
        if(n_type_used>1)
           SID_trap_error("An invalid number of particle types (%d) are being used in make_map.",ERROR_LOGIC,n_type_used);
        mq->x=(GBPREAL *)ADaPS_fetch(plist->data,"x_%s",plist->species[i_type]);
        mq->y=(GBPREAL *)ADaPS_fetch(plist->data,"y_%s",plist->species[i_type]);
        mq->z=(GBPREAL *)ADaPS_fetch(plist->data,"z_%s",plist->species[i_type]);
        if(ADaPS_exist(plist->data,"r_smooth_%s",plist->species[i_type]))
           mq->h_smooth=(float *)ADaPS_fetch(plist->data,"r_smooth_%s",plist->species[i_type]);
        /*
        else if(ADaPS_exist(plist->data,"rho_%s",plist->species[i_type])){
          (*h_smooth)=NULL;
          if((*rho)==NULL)
            (*rho)=(float *)ADaPS_fetch(plist->data,"rho_%s",plist->species[i_type]);
          if(ADaPS_exist(plist->data,"M_%s",plist->species[i_type])){
            if((*mass)==NULL)
               (*mass)=(double *)ADaPS_fetch(plist->data,"M_%s",plist->species[i_type]);
          }
          else{
            (*h_smooth)=NULL;
            if(ADaPS_exist(plist->data,"mass_array_%s",plist->species[i_type])){
              if((*mass_array)==NULL)
                 (*mass_array)=((double *)ADaPS_fetch(plist->data,"M_%s",plist->species[i_type]))[0];
              else
                mass_array=1.;
          }
        }
        */
        else
           SID_trap_error("No smoothing lengths available for type={%s} in make_map.",ERROR_LOGIC,plist->species[i_type]);
     }
  }
  SID_log("%zd eligible for rendering...Done.",SID_LOG_CLOSE,mq->n_particles);

}

void set_particle_map_quantities(map_quantities_info *mq,int mode,size_t i_particle,
                                 float *x_i,
                                 float *y_i,
                                 float *z_i,
                                 float *h_i,
                                 float *v_i,
                                 float *w_i);
void set_particle_map_quantities(map_quantities_info *mq,int mode,size_t i_particle,
                                 float *x_i,
                                 float *y_i,
                                 float *z_i,
                                 float *h_i,
                                 float *v_i,
                                 float *w_i){
  
  (*x_i)=mq->x[i_particle];
  (*y_i)=mq->y[i_particle];
  (*z_i)=mq->z[i_particle];
  if(mode==TRUE){
     (*h_i)=mq->h_smooth[i_particle];
     // Set the particle weighting
     switch(mq->w_mode){
        case MAKE_MAP_MODE_RHO:
           (*w_i)=mq->rho[i_particle];
           if(mq->transfer_rho!=NULL){
              switch(mq->flag_transfer_rho_log){
                 case TRUE:
                    (*w_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,take_log10((double)(*w_i)))));
                    break;
                 case FALSE:
                    (*w_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,(double)(*w_i))));
                    break;
              }
              if(!mq->flag_comoving)
                 (*w_i)*=mq->inv_expansion_factor_cubed;
           }
           break;
        case MAKE_MAP_NO_WEIGHTING:
           break;
        default:
           SID_trap_error("Unknown w_mode (%d) in make_map.",ERROR_LOGIC,mq->v_mode);
           break;
     }
     // Set the particle value
     switch(mq->v_mode){
        case MAKE_MAP_MODE_RHO:
           (*v_i)=mq->rho[i_particle];
           if(mq->transfer_rho!=NULL){
              switch(mq->flag_transfer_rho_log){
                 case TRUE:
                    (*v_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,take_log10((double)(*v_i)))));
                    break;
                 case FALSE:
                    (*v_i)*=MAX(0.,MIN(1.,interpolate(mq->transfer_rho,(double)(*v_i))));
                    break;
              }
              if(!mq->flag_comoving)
                 (*v_i)*=mq->inv_expansion_factor_cubed;
           }
           break;
        case MAKE_MAP_MODE_SIGMA:
           (*v_i)=mq->sigma[i_particle];
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

float compute_f_stretch(double d_xyz,float z_i,int flag_plane_parallel);
float compute_f_stretch(double d_xyz,float z_i,int flag_plane_parallel){
   float f_i;
   switch(flag_plane_parallel){
      case FALSE:
         if(z_i>0.)
            f_i=(float)d_xyz/(float)z_i;
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
//   needed to place the camera at (0,0,-d_xyz)
//   with the object at (0,0,0)
void compute_perspective_transformation(double  x_o,
                                        double  y_o,
                                        double  z_o,
                                        double  x_c,
                                        double  y_c,
                                        double  z_c,
                                        double  stereo_offset,
                                        double *d_xy,
                                        double *d_xyz,
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
                                        double  stereo_offset,
                                        double *d_xy,
                                        double *d_xyz,
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
  d_x_o   = x_o-x_c;
  d_y_o   = y_o-y_c;
  d_z_o   = z_o-z_c;
  (*d_xy) = sqrt(pow(d_x_o,2.)+pow(d_y_o,2.));
  (*d_xyz)= sqrt(pow(d_x_o,2.)+pow(d_y_o,2.)+pow(d_z_o,2.));
  (*x_hat)= d_y_o/(*d_xy);
  (*y_hat)=-d_x_o/(*d_xy);
  (*z_hat)= 0.;
  (*theta)= acos(d_z_o/(*d_xyz));
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
     SID_log("Forcing theta=0 in stereo projection.",SID_LOG_COMMENT);
     Dx_stereo=stereo_offset*cos(theta_roll_stereo)*d_y_o/(*d_xy);
     Dy_stereo=stereo_offset*cos(theta_roll_stereo)*d_x_o/(*d_xy);
     if(theta_roll_stereo>0.)
        Dz_stereo=sqrt(stereo_offset*stereo_offset-Dx_stereo*Dx_stereo-Dy_stereo*Dy_stereo);
     else if(theta_roll_stereo==0.)
        Dz_stereo=0.;
     else
        Dz_stereo=-sqrt(stereo_offset*stereo_offset-Dx_stereo*Dx_stereo-Dy_stereo*Dy_stereo);
     (*x_o_out)+=Dx_stereo;
     (*y_o_out)+=Dy_stereo;
     (*z_o_out)+=Dz_stereo;
     (*x_c_out)+=Dx_stereo;
     (*y_c_out)+=Dy_stereo;
     (*z_c_out)+=Dz_stereo;
  }
}

void init_make_map(plist_info  *plist,
                   char        *parameter,
                   ADaPS       *transfer_list,
                   double       x_o,
                   double       y_o,
                   double       z_o,
                   double       x_c,
                   double       y_c,
                   double       z_c,
                   double       box_size,
                   double       FOV_x,
                   double       FOV_y,
                   double       xmin,
                   double       ymin,
                   double       pixel_size_x,
                   double       pixel_size_y,
                   double       radius_kernel_norm,
                   int          nx,
                   int          ny,
                   int          mode,
                   double       expansion_factor,
                   double       near_field,
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
                   size_t      **z_index,
                   int          *i_x_min_local_return,
                   int          *i_x_max_local_return,
                   size_t       *n_particles);
void init_make_map(plist_info  *plist,
                   char        *parameter,
                   ADaPS       *transfer_list,
                   double       x_o,
                   double       y_o,
                   double       z_o,
                   double       x_c,
                   double       y_c,
                   double       z_c,
                   double       box_size,
                   double       FOV_x,
                   double       FOV_y,
                   double       xmin,
                   double       ymin,
                   double       pixel_size_x,
                   double       pixel_size_y,
                   double       radius_kernel_norm,
                   int          nx,
                   int          ny,
                   int          mode,
                   double       expansion_factor,
                   double       near_field,
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
  double   d_xy;
  double   d_xyz;
  double   d_hat;
  double   particle_radius;
  double   x_tmp,y_tmp,z_tmp;
  interp_info *transfer;
  double       transfer_val;
  int          flag_log;
  double       z_test;
  int          flag_use_Gadget;
  
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
  init_particle_map_quantities(&mq,plist,transfer_list,parameter,flag_comoving,expansion_factor);
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
  compute_perspective_transformation(x_o_in,
                                     y_o_in,
                                     z_o_in,
                                     x_c_in,
                                     y_c_in,
                                     z_c_in,
                                     stereo_offset,
                                     &d_xy,
                                     &d_xyz,
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

  // Determine how many particles are contributing to each 
  //    column of the image
  float   x_i;
  float   y_i;
  float   z_i;
  float   h_i;
  float   f_i;
  float   v_i;
  float   w_i;
  int     i_x;
  size_t *n_column;
  size_t  n_visible;
  size_t  n_visible_local=0;
  n_column=(size_t *)SID_calloc(sizeof(size_t)*nx);
  SID_log("Compute slab domain decomposition...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
    if(ptype_used[i_type] && ADaPS_exist(plist->data,"n_%s",plist->species[i_type])){
      n_particles_species=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_type]))[0];
      for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){

         // Set the preoperties of the particle to be mapped
         set_particle_map_quantities(&mq,TRUE,k_particle,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

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
                            d_xyz,
                            stereo_offset,
                            theta,
                            theta_roll,
                            box_size,
                            expansion_factor,
                            flag_comoving,
                            flag_force_periodic);


         if(z_i>0){
            double radius2_norm;
            double radius_kernel;
            double part_pos_x;
            double part_pos_y;
            int    kx_min;
            int    kx_max;
            int    ky_min;
            int    ky_max;
            int    n_ky;
            f_i=(float)compute_f_stretch(d_xyz,z_i,flag_plane_parallel);
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
                            radius_kernel_norm,
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
            n_ky  =ky_max-ky_min+1;
            for(i_x=kx_min;i_x<=kx_max;i_x++)
               n_column[i_x]+=n_ky;
            if(kx_min<=kx_max)
               n_visible_local++;
         }
      }
      j_particle+=n_particles_species;
    }
  }
  SID_Allreduce(SID_IN_PLACE,    n_column,  nx,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  SID_Allreduce(&n_visible_local,&n_visible,1, SID_SIZE_T,SID_SUM,SID.COMM_WORLD);

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
     for(;i_x<nx-1 && n_column[i_x]<n_target;) i_x++;
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
    if(ptype_used[i_type] && ADaPS_exist(plist->data,"n_%s",plist->species[i_type])){
      n_particles_species=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_type]))[0];
      for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){

         // Set the preoperties of the particle to be mapped
         set_particle_map_quantities(&mq,TRUE,k_particle,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

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
                            d_xyz,
                            stereo_offset,
                            theta,
                            theta_roll,
                            box_size,
                            expansion_factor,
                            flag_comoving,
                            flag_force_periodic);

         if(z_i>0){
            double radius2_norm;
            double radius_kernel;
            double part_pos_x;
            double part_pos_y;
            int    kx_min;
            int    kx_max;
            int    ky_min;
            int    ky_max;
            f_i=(float)compute_f_stretch(d_xyz,z_i,flag_plane_parallel);
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
                            radius_kernel_norm,
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
            for(i_rank=0;i_rank<SID.n_proc;i_rank++){
               if(kx_min<=i_x_max_rank[i_rank] && kx_max>=i_x_min_rank[i_rank])
                  n_rank_local[i_rank]++;
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
  size_t i_particle_rank;
  size_t j_particle_rank=0;
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
     int rank_to;
     int rank_from;
     set_exchange_ring_ranks(&rank_to,&rank_from,i_rank);
     i_particle_rank=0;
     for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
       if(ptype_used[i_type]){
         n_particles_species=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_type]))[0];
         for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){

            // Set the preoperties of the particle to be mapped
            set_particle_map_quantities(&mq,TRUE,k_particle,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

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
                               d_xyz,
                               stereo_offset,
                               theta,
                               theta_roll,
                               box_size,
                               expansion_factor,
                               flag_comoving,
                               flag_force_periodic);

            if(z_i>0){
               double radius2_norm;
               double radius_kernel;
               double part_pos_x;
               double part_pos_y;
               int    kx_min;
               int    kx_max;
               int    ky_min;
               int    ky_max;
               f_i=(float)compute_f_stretch(d_xyz,z_i,flag_plane_parallel);
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
                               radius_kernel_norm,
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
               if(kx_min<=i_x_max_rank[rank_to] && kx_max>=i_x_min_rank[rank_to]){
                  x_buffer[i_particle_rank]=x_i;
                  y_buffer[i_particle_rank]=y_i;
                  z_buffer[i_particle_rank]=z_i;
                  h_buffer[i_particle_rank]=h_i;
                  f_buffer[i_particle_rank]=f_i;
                  w_buffer[i_particle_rank]=w_i;
                  v_buffer[i_particle_rank]=v_i;
                  i_particle_rank++;
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
  free_particle_map_quantities(&mq);

  (*i_x_min_local_return)=i_x_min_local;
  (*i_x_max_local_return)=i_x_max_local;

  SID_log("Done.",SID_LOG_CLOSE);
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
  float     *weight;
  float     *value;
  size_t    *z_index;
  int        kx_min;
  int        kx_max;
  int        ky_min;
  int        ky_max;
  int        kz_min;
  int        kz_max;
  double    *numerator;
  double    *denominator;
  double    *z_image;
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
  double     min_image;
  double     max_image;
  double     min_z_image;
  double     max_z_image;
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
  double     radius_kernel_norm;
  double     radius_kernel_norm2;

  double     d_o,d_x_o,d_y_o,d_z_o;

  double       x_o;
  double       y_o;
  double       z_o;
  double       x_c;
  double       y_c;
  double       z_c;
  double       FOV_x;
  double       FOV_y;
  double       box_size;
  int          nx;
  int          ny;
  double      *image;
  char        *mask;
  int          mode;
  char        *parameter;
  plist_info  *plist;
  int          flag_comoving;
  int          flag_force_periodic;
  double       expansion_factor;
  ADaPS       *transfer;
  double       h_Hubble;
  double       near_field;
  double       stereo_offset;
  int          i_x,i_y;
  int         *mask_buffer;

  int          i_image;
  int          camera_mode;

  // Create an absorption look-up table to speed things up
  int          flag_add_absorption;
  double       absorption_coefficient;
  double       inv_kappa_absorption;
  double      *column_depth;
  double      *x_abs;
  double      *y_abs;
  size_t       n_abs=400;
  int          i_abs;
  double       tau_max=10.;
  interp_info *abs_interp;
  double       rho_abs_lo;
  double       rho_abs_hi;
  double       inv_kappa_abs_lo;
  double       inv_kappa_abs_hi;
  interp_info *inv_kappa_interp;

  if(render->flag_add_absorption){
     flag_add_absorption       =TRUE;
     if(render->kappa_transfer){
        int     n_temp;
        int     i_temp;
        double *inv_y_temp;
        n_temp    =render->kappa_transfer->n;
        inv_y_temp=(double *)SID_malloc(sizeof(double)*n_temp);
        for(i_temp=0;i_temp<n_temp;i_temp++)
           inv_y_temp[i_temp]=-render->kappa_transfer->y[i_temp];
        for(i_temp=0;i_temp<n_temp;i_temp++)
           inv_y_temp[i_temp]=1./render->kappa_transfer->y[i_temp];
        rho_abs_lo      =take_alog10(render->kappa_transfer->x[0]);
        rho_abs_hi      =take_alog10(render->kappa_transfer->x[n_temp-1]);
        inv_kappa_abs_lo=take_alog10(inv_y_temp[0]);
        inv_kappa_abs_hi=take_alog10(inv_y_temp[n_temp-1]);
        rho_abs_lo      =render->kappa_transfer->x[0];
        rho_abs_hi      =render->kappa_transfer->x[n_temp-1];
        inv_kappa_abs_lo=inv_y_temp[0];
        inv_kappa_abs_hi=inv_y_temp[n_temp-1];
        init_interpolate(render->kappa_transfer->x,inv_y_temp,(size_t)n_temp,gsl_interp_linear,&inv_kappa_interp);
        /*
        if(SID.I_am_Master){
          double x_test;
          for(x_test=rho_abs_lo;x_test<rho_abs_hi;x_test+=(rho_abs_hi-rho_abs_lo)/100.)
            fprintf(stderr,"%le %le %le\n",x_test,interpolate(render->kappa_transfer,x_test),interpolate(inv_kappa_interp,x_test));
        }
        */
        SID_free(SID_FARG inv_y_temp);
     }
     else{
        inv_kappa_interp      =NULL;
        absorption_coefficient=render->kappa_absorption;
        inv_kappa_absorption  =1./absorption_coefficient;
     }
     x_abs   =(double *)SID_malloc(sizeof(double)*n_abs);
     y_abs   =(double *)SID_malloc(sizeof(double)*n_abs);
     x_abs[0]=0.;
     for(i_abs=1;i_abs<(n_abs-1);i_abs++)
        x_abs[i_abs]=(float)i_abs*((float)tau_max/(float)(n_abs-1));
     x_abs[n_abs-1]=tau_max;
     for(i_abs=0;i_abs<n_abs;i_abs++)
        y_abs[i_abs]=exp(-x_abs[i_abs]);
     init_interpolate(x_abs,y_abs,n_abs,gsl_interp_linear,&abs_interp);
     SID_free(SID_FARG x_abs);
     SID_free(SID_FARG y_abs);
  }
  else{
     flag_add_absorption   =FALSE;
     absorption_coefficient=0.;
     inv_kappa_absorption  =0.;
  }

  plist   =&(render->plist);
  x_o     =render->camera->perspective->p_o[0];
  y_o     =render->camera->perspective->p_o[1];
  z_o     =render->camera->perspective->p_o[2];
  d_o     =render->camera->perspective->d_o;
  x_c     =render->camera->perspective->p_c[0];
  y_c     =render->camera->perspective->p_c[1];
  z_c     =render->camera->perspective->p_c[2];
  nx      =render->camera->width;
  ny      =render->camera->height;
  mode    =render->mode;
  camera_mode        =render->camera->camera_mode;
  flag_comoving      =render->flag_comoving;
  flag_force_periodic=render->flag_force_periodic;
  expansion_factor   =render->camera->perspective->time;
  box_size        =((double *)ADaPS_fetch(plist->data,"box_size"))[0];
  h_Hubble        =render->h_Hubble;
  near_field      =render->near_field*d_o;
  SID_log("near field =%le [Mpc/h]",SID_LOG_COMMENT,near_field*h_Hubble/M_PER_MPC);
  if(nx>=ny){
    FOV_y=render->camera->perspective->FOV;
    FOV_x=FOV_y*(double)nx/(double)ny;    
  }
  else{
    FOV_x=render->camera->perspective->FOV;
    FOV_y=FOV_x*(double)nx/(double)ny;        
  }

  if(check_mode_for_flag(camera_mode,CAMERA_STEREO)){
    n_images=6;
    i_image =2; // Skip mono images
  }
  else{
    n_images=2;
    i_image =0;
  }

  for(;i_image<n_images;i_image++){

    switch(i_image){
      // Mono RGB
      case 0:
        parameter=render->camera->RGB_param;
        transfer =render->camera->RGB_transfer;
        image    =render->camera->image_RGB->values;
        z_image  =NULL;
        mask     =render->camera->mask_RGB;
        if(nx>=ny){
          FOV_y=render->camera->perspective->FOV;
          FOV_x=FOV_y*(double)nx/(double)ny;
        }
        else{
          FOV_x=render->camera->perspective->FOV;
          FOV_y=FOV_x*(double)nx/(double)ny;
        }
        stereo_offset=0.;
        break;
      // Mono luminance
      case 1:
        parameter=render->camera->Y_param;
        transfer =render->camera->Y_transfer;
        image    =render->camera->image_Y->values;
        if(render->camera->image_Z!=NULL)
           z_image=render->camera->image_Z->values;
        else
           z_image=NULL;
        mask=render->camera->mask_Y;
        break;
      // Left RGB
      case 2:
        parameter     =render->camera->RGB_param;
        transfer      =render->camera->RGB_transfer;
        image         =render->camera->image_RGB_left->values;
        z_image       =NULL;
        mask          =render->camera->mask_RGB_left;
        stereo_offset =-d_o/render->camera->stereo_ratio;
        break;
      // Left luminance
      case 3:
        parameter=render->camera->Y_param;
        transfer =render->camera->Y_transfer;
        image    =render->camera->image_Y_left->values;
        if(render->camera->image_Z_left!=NULL)
           z_image=render->camera->image_Z_left->values;
        else
           z_image=NULL;
        mask=render->camera->mask_Y_left;
        break;
      // Right RGB
      case 4:
        parameter     =render->camera->RGB_param;
        transfer      =render->camera->RGB_transfer;
        image         =render->camera->image_RGB_right->values;
        z_image       =NULL;
        mask          =render->camera->mask_RGB_right;
        stereo_offset =d_o/render->camera->stereo_ratio;
        break;
      // Right luminance
      case 5:
        parameter=render->camera->Y_param;
        transfer =render->camera->Y_transfer;
        image    =render->camera->image_Y_right->values;
        if(render->camera->image_Y_right!=NULL)
           z_image=render->camera->image_Z_right->values;
        else
           z_image=NULL;
        mask=render->camera->mask_Y_right;
        break;
    }

    SID_log("Projecting {%s} to a %dx%d pixel array...",SID_LOG_OPEN|SID_LOG_TIMER,parameter,nx,ny);

    // Set physical image-plane domain
    xmin  = -FOV_x/2.; // Things will be centred on (x_o,y_o,z_o) later
    ymin  = -FOV_y/2.; // Things will be centred on (x_o,y_o,z_o) later

    // Apply frustrum correction
    if(stereo_offset<0.)
       xmin-=stereo_offset;
  
    // Compute image scales
    n_pixels    =nx*ny;
    pixel_size_x=FOV_x/(double)nx;
    pixel_size_y=FOV_y/(double)ny;
    pixel_size  =0.5*(pixel_size_x+pixel_size_y);
    pixel_area  =pixel_size_x*pixel_size_y;
    if(fabs((pixel_size_x-pixel_size_y)/pixel_size_x)>1e-4)
      SID_log_warning("pixels are not square by %7.3f%%",0,fabs((pixel_size_x-pixel_size_y)/pixel_size_x)*1e2);

    // Generate the smoothing kernal
    set_sph_kernel(plist,SPH_KERNEL_GADGET|SPH_KERNEL_2D);
    kernel_radius      =(double *)ADaPS_fetch(plist->data, "sph_kernel_radius");
    kernel_table       =(double *)ADaPS_fetch(plist->data, "sph_kernel_2d");  
    kernel_table_avg   =((double *)ADaPS_fetch(plist->data,"sph_kernel_2d_dA"))[0];  
    radius_kernel_norm =kernel_radius[N_KERNEL_TABLE];
    radius_kernel_norm2=radius_kernel_norm*radius_kernel_norm;

    double x_o_in=x_o;
    double y_o_in=y_o;
    double z_o_in=z_o;
    double x_c_in=x_c;
    double y_c_in=y_c;
    double z_c_in=z_c;
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
    double d_xy;
    double d_xyz;
    compute_perspective_transformation(x_o_in,
                                       y_o_in,
                                       z_o_in,
                                       x_c_in,
                                       y_c_in,
                                       z_c_in,
                                       stereo_offset,
                                       &d_xy,
                                       &d_xyz,
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
    init_make_map(plist,
                  parameter,
                  transfer,
                  x_o,y_o,z_o,
                  x_c,y_c,z_c,
                  box_size,FOV_x,FOV_y,
                  xmin,ymin,
                  pixel_size_x,pixel_size_y,
                  radius_kernel_norm,
                  nx,ny,
                  mode,
                  expansion_factor,
                  near_field,
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
                  &z_index,
                  &i_x_min_local,
                  &i_x_max_local,
                  &n_particles);

    // Allocate and initialize image arrays
    if(flag_add_absorption)
       column_depth=(double *)SID_malloc(sizeof(double)*n_pixels);
    else
       column_depth=NULL;
    if(flag_weigh)
       denominator=(double *)SID_malloc(sizeof(double)*n_pixels);
    else
       denominator=NULL;
    numerator=image;
    for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
      image[i_pixel]    =0.;
      numerator[i_pixel]=0.;
      if(column_depth!=NULL)
        column_depth[i_pixel]=0.;
      if(denominator!=NULL)
        denominator[i_pixel]=0.;
      if(z_image!=NULL)
        z_image[i_pixel] =0.;
      mask[i_pixel]=FALSE;
    }

    d_x_o=x_o-x_c;
    d_y_o=y_o-y_c;
    d_z_o=z_o-z_c;
    d_o  =sqrt(pow(d_x_o,2.)+
               pow(d_y_o,2.)+
               pow(d_z_o,2.));

    // Report absorption statistics
    if(render->kappa_transfer!=NULL){
       float rho_min,rho_mean,rho_max;
       if(denominator!=NULL){
          calc_min_global( weight,&rho_min, n_particles,SID_FLOAT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
          calc_max_global( weight,&rho_max, n_particles,SID_FLOAT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
          calc_mean_global(weight,&rho_mean,n_particles,SID_FLOAT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
       }
       else{
          calc_min_global( value,&rho_min, n_particles,SID_FLOAT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
          calc_max_global( value,&rho_max, n_particles,SID_FLOAT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
          calc_mean_global(value,&rho_mean,n_particles,SID_FLOAT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
       }
       SID_log("Absorption statistics:",SID_LOG_OPEN);
       SID_log("rho_min    =%le",     SID_LOG_COMMENT,rho_min);
       SID_log("rho_max    =%le",     SID_LOG_COMMENT,rho_max);
       SID_log("rho_mean   =%le",     SID_LOG_COMMENT,rho_mean);
       SID_log("table range=%le->%le",SID_LOG_COMMENT,rho_abs_lo,rho_abs_hi);
       SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
    }  

    // Perform projection
    size_t        ii_particle;
    pcounter_info pcounter;
    SID_init_pcounter(&pcounter,n_particles,10);
    SID_log("Performing projection...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(ii_particle=0;ii_particle<n_particles;ii_particle++){
      i_particle=z_index[ii_particle];
      z_i       =(double)z[i_particle];

      part_h_z=(double)h_smooth[i_particle];

      if(z_i>near_field){

        // Set pixel space ranges and positions
        /*
        set_pixel_space(h_smooth[i_particle],
                        x[i_particle],
                        y[i_particle],
                        f_stretch[i_particle],
                        xmin,
                        ymin,
                        FOV_x,
                        FOV_y,
                        pixel_size_x,
                        pixel_size_y,
                        radius_kernel_norm,
                        &radius2_norm,
                        &radius_kernel,
                        &part_pos_x,
                        &part_pos_y,
                        &kx_min,
                        &kx_max,
                        &ky_min,
                        &ky_max);
        */

        part_h_xy    =part_h_z*f_stretch[i_particle];
        radius2_norm =1./(part_h_xy*part_h_xy);
        part_pos_x   =(double)(x[i_particle]*f_stretch[i_particle]);
        part_pos_y   =(double)(y[i_particle]*f_stretch[i_particle]);
        radius_kernel=radius_kernel_norm*part_h_xy;
        kx_min=(int)((part_pos_x-radius_kernel-xmin)/pixel_size_x);
        kx_max=(int)((part_pos_x+radius_kernel-xmin)/pixel_size_x+ONE_HALF);
        ky_min=(int)((part_pos_y-radius_kernel-ymin)/pixel_size_y);
        ky_max=(int)((part_pos_y+radius_kernel-ymin)/pixel_size_y+ONE_HALF);

        // Set the particle values and weights
        double w_i;
        double v_i;
        double vw_i;
        v_i=(double)value[i_particle];
        if(denominator!=NULL){
           w_i =(double)weight[i_particle];
           vw_i=v_i*w_i;
        }

        // Loop over the kernal
        for(kx=kx_min,pixel_pos_x=xmin+(kx_min+0.5)*pixel_size_x;kx<=kx_max;kx++,pixel_pos_x+=pixel_size_x){
          if(kx>=0 && kx<nx){
            for(ky=ky_min,pixel_pos_y=ymin+(ky_min+0.5)*pixel_size_y;ky<=ky_max;ky++,pixel_pos_y+=pixel_size_y){
              if(ky>=0 && ky<ny){
                radius2=
                  (pixel_pos_x-part_pos_x)*(pixel_pos_x-part_pos_x)+
                  (pixel_pos_y-part_pos_y)*(pixel_pos_y-part_pos_y);
                radius2*=radius2_norm;
                // Construct image here
                if(radius2<radius_kernel_norm2){
                  pos    =ky+kx*ny;
                  f_table=sqrt(radius2);
                  i_table=(int)(f_table*(double)N_KERNEL_TABLE);
                  kernel =kernel_table[i_table]+
                    (kernel_table[i_table+1]-kernel_table[i_table])*
                    (f_table-kernel_radius[i_table])*(double)N_KERNEL_TABLE;
                  if(denominator!=NULL){
                    switch(flag_add_absorption){
                       case TRUE:
                          double absorption;
                          double tau;
                          double dtau;
                          if(inv_kappa_interp==NULL)
                             dtau=w_i*kernel*inv_kappa_absorption;
                          else if(w_i<rho_abs_lo)
                             dtau=w_i*kernel*inv_kappa_abs_lo;
                          else if(w_i>rho_abs_hi)
                             dtau=w_i*kernel*inv_kappa_abs_hi;
                          else
                             dtau=w_i*kernel*interpolate(inv_kappa_interp,w_i);
                          tau=(double)column_depth[pos];
                          if(tau>tau_max)
                             absorption=0.;
                          else if(tau<=0.)
                             absorption=1.;
                          else
                             absorption=interpolate(abs_interp,tau);
                          //fprintf(stderr,"%e %le %le %le\n",w_i,interpolate(inv_kappa_interp,w_i),dtau,absorption);
                          column_depth[pos]+=dtau;
                          numerator[pos]   +=vw_i*kernel*absorption;
                          denominator[pos] += w_i*kernel*absorption;
                          if(z_image!=NULL)
                             z_image[pos]+=z_i*w_i*kernel*absorption;
                          break;
                       case FALSE:
                          numerator[pos]  +=vw_i*kernel;
                          denominator[pos]+=w_i*kernel;
                          if(z_image!=NULL)
                             z_image[pos]+=z_i*w_i*kernel;
                          break;
                    }
                  }
                  else{
                    switch(flag_add_absorption){
                       case TRUE:
                          double absorption;
                          double tau;
                          double dtau;
                          if(inv_kappa_interp==NULL)
                             dtau=v_i*kernel*inv_kappa_absorption;
                          else if(v_i<rho_abs_lo)
                             dtau=v_i*kernel*inv_kappa_abs_lo;
                          else if(v_i>rho_abs_hi)
                             dtau=v_i*kernel*inv_kappa_abs_hi;
                          else
                             dtau=v_i*kernel*interpolate(inv_kappa_interp,v_i);
                          tau=(double)column_depth[pos];
                          if(tau>tau_max)
                             absorption=0.;
                          else if(tau<=0.)
                             absorption=1.;
                          else 
                             absorption=interpolate(abs_interp,tau);
                          column_depth[pos]+=dtau;
                          numerator[pos]   +=v_i*kernel*absorption;
                          if(z_image!=NULL)
                             z_image[pos]+=z_i*v_i*kernel*absorption;
                          break;
                       case FALSE:
                          numerator[pos]+=v_i*kernel;
                          if(z_image!=NULL)
                             z_image[pos]+=z_i*v_i*kernel;
                          break;
                    }
                  }
                  mask[pos] =TRUE;
                }
              }
            }
          }
        }
      }
      SID_check_pcounter(&pcounter,ii_particle);
    }
    SID_Barrier(SID.COMM_WORLD);
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Image normalization, etc...",SID_LOG_OPEN|SID_LOG_TIMER);

    // Add results from all ranks if this is being run in parallel
    //   First, clear the parts of the image not in this rank's domain ...
    for(kx=0;kx<nx;kx++){
       if(!(kx>=i_x_min_local && kx<=i_x_max_local)){
          for(ky=0;ky<ny;ky++){
             pos=ky+kx*ny;
             numerator[pos]  =0.;
             if(denominator!=NULL)
                denominator[pos]=0.;
             if(z_image!=NULL)
                z_image[pos]=0.;
             if(column_depth!=NULL)
                column_depth[pos]=0.;
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
    SID_Allreduce(SID_IN_PLACE,numerator,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
    if(denominator!=NULL)
      SID_Allreduce(SID_IN_PLACE,denominator,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);
    if(z_image!=NULL)
      SID_Allreduce(SID_IN_PLACE,z_image,n_pixels,SID_DOUBLE,SID_SUM,SID.COMM_WORLD);

    // Normalize image (if needed)
    if(denominator!=NULL){
      for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
        if(mask[i_pixel])
          image[i_pixel]/=denominator[i_pixel];
      }
    }

    // Normalize z-frame (if needed)
    if(z_image!=NULL){
      if(denominator!=NULL){
        for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
            z_image[i_pixel]=z_image[i_pixel]/denominator[i_pixel];
          else
            z_image[i_pixel]=0.;
        }
      }
      else{
        for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
            z_image[i_pixel]=z_image[i_pixel]/numerator[i_pixel];
          else
            z_image[i_pixel]=0.;
        }
      }
    }

    // Take log_10 (if needed)
    if(check_mode_for_flag(mode,MAKE_MAP_LOG)){
      for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
        if(mask[i_pixel])
          image[i_pixel]=take_log10(image[i_pixel]);
        else
          image[i_pixel]=LOG_ZERO;
      }
    }

    // Compute some image statistics
    for(i_pixel=0,n_unmasked=0;i_pixel<n_pixels;i_pixel++) 
      if(mask[i_pixel]) 
        n_unmasked++;
    calc_min(image,&min_image, n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    calc_max(image,&max_image, n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    if(z_image!=NULL){
      calc_min(z_image,&min_z_image, n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
      calc_max(z_image,&max_z_image, n_pixels,SID_DOUBLE,CALC_MODE_DEFAULT);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  
    // Report image statistics
    if(n_unmasked>0){
      SID_log("Image statistics:   (min,max,coverage)=(%8.3le,%8.3le,%3d%%)",SID_LOG_COMMENT,min_image,max_image,(int)(100.*n_unmasked/n_pixels));
      if(z_image!=NULL)
        SID_log("Z-frame statistics: (min,max)         =(%8.3le,%8.3le) [Mpc/h]",SID_LOG_COMMENT,h_Hubble*min_z_image/M_PER_MPC,h_Hubble*max_z_image/M_PER_MPC);
    }
    else
      SID_out("Image is empty.",SID_LOG_COMMENT);

    // Clean-up
    SID_free(SID_FARG x);
    SID_free(SID_FARG y);
    SID_free(SID_FARG z);
    SID_free(SID_FARG h_smooth);
    SID_free(SID_FARG f_stretch);
    SID_free(SID_FARG value);
    SID_free(SID_FARG z_index);
    SID_free(SID_FARG weight);
    if(denominator!=NULL)
      SID_free(SID_FARG denominator);
    if(column_depth!=NULL)
      SID_free(SID_FARG column_depth);
  
    SID_log("Done.",SID_LOG_CLOSE);
  }

  if(flag_add_absorption){
     free_interpolate(SID_FARG abs_interp);
     if(inv_kappa_interp!=NULL)
        free_interpolate(SID_FARG inv_kappa_interp);
  }

}

