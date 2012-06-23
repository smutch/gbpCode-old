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

   // Apply stereo offset
   (*x_i)+=2.*stereo_offset;
   
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
   double        mass_array;
   float        *rho;
   interp_info  *transfer_rho;
   int           flag_transfer_rho_log;
   float        *sigma;
   interp_info  *transfer_sigma;
   int           flag_transfer_sigma_log;
   int           flag_comoving;
   double        inv_expansion_factor_cubed;
   int           v_mode;
   int           w_mode;
};

void init_particle_map_quantities(map_quantities_info *mq,plist_info *plist,ADaPS *transfer_list,char *parameter,int flag_comoving,double expansion_factor);
void init_particle_map_quantities(map_quantities_info *mq,plist_info *plist,ADaPS *transfer_list,char *parameter,int flag_comoving,double expansion_factor){

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
                   double       FOV,
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
                   size_t      *n_particles);
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
                   double       FOV,
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
                   size_t      *n_particles){
  int      ptype_used[N_GADGET_TYPE];
  int      i_type;
  float   *rho;
  double  *mass;
  float   *sigma_v;
  size_t   i_particle;
  size_t   j_particle;
  size_t   k_particle;
  double   mass_array;
  size_t   n_particles_species;
  float    *x_temp;
  float    *y_temp;
  float    *z_temp;
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
  double   half_box_size;
  double   particle_radius;
  double   x_tmp,y_tmp,z_tmp;
  interp_info *transfer;
  double       transfer_val;
  int          flag_log;
  double       z_test;
  int          flag_use_Gadget;
  
  SID_log("Initializing projection-space...",SID_LOG_OPEN|SID_LOG_TIMER);

  half_box_size=0.5*box_size;

  // Plane parallel projection?
  int flag_plane_parallel;
  if(check_mode_for_flag(camera_mode,CAMERA_PLANE_PARALLEL))
     flag_plane_parallel=TRUE;
  else
     flag_plane_parallel=FALSE;

  // Initialize the mapping quantities
  map_quantities_info mq;
  init_particle_map_quantities(&mq,plist,transfer_list,parameter,flag_comoving,expansion_factor);

  // Compute the angle and the axis of rotation
  //   needed to place the camera at (0,0,-d_o)
  //   with the object at (0,0,0)
  x_o-=stereo_offset;
  y_o-=stereo_offset;
  z_o-=stereo_offset;
  x_c-=stereo_offset;
  y_c-=stereo_offset;
  z_c-=stereo_offset;
  d_x_o= x_o-x_c;
  d_y_o= y_o-y_c;
  d_z_o= z_o-z_c;
  d_o  = sqrt(pow(d_x_o,2.)+pow(d_y_o,2.));
  x_hat= d_y_o/d_o;
  y_hat=-d_x_o/d_o;
  z_hat= 0.;
  d_o  = sqrt(pow(d_x_o,2.)+pow(d_y_o,2.)+pow(d_z_o,2.));
  theta= acos(d_z_o/d_o);
  if(sqrt(d_x_o*d_x_o+d_y_o*d_y_o)>0.){
    theta_roll=acos(-d_y_o/sqrt(d_x_o*d_x_o+d_y_o*d_y_o));
    if(d_x_o<0.)
      theta_roll=TWO_PI-theta_roll;
  }
  else
    theta_roll=0.;
  
  // Initialize the local minima of each array we will exchange later
  float x_min=FLT_MAX;
  float y_min=FLT_MAX;
  float z_min=FLT_MAX;
  float h_min=FLT_MAX;
  float f_min=FLT_MAX;
  float v_min=FLT_MAX;
  float w_min=FLT_MAX;

  // First, compute the positions along the line of sight
  //    so we can perform the domain decomposition
  float  x_i;
  float  y_i;
  float  z_i;
  float  h_i;
  float  f_i;
  float  v_i;
  float  w_i;
  float *z_decomp;
  int   *keep;
  int    flag_init=TRUE;
  size_t n_visible_local=0;
  size_t n_visible=0;
  int    i_rank;
  z_decomp=(float *)SID_malloc(sizeof(float)*(*n_particles));
  for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
    if(ptype_used[i_type] && ADaPS_exist(plist->data,"n_%s",plist->species[i_type])){
      n_particles_species=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_type]))[0];
      // Fetch coordinates
      x_temp=(float *)ADaPS_fetch(plist->data,"x_%s",plist->species[i_type]);
      y_temp=(float *)ADaPS_fetch(plist->data,"y_%s",plist->species[i_type]);
      z_temp=(float *)ADaPS_fetch(plist->data,"z_%s",plist->species[i_type]);
      for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){

         // Set the preoperties of the particle to be mapped
         set_particle_map_quantities(&mq,FALSE,i_particle,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);

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
                            flag_comoving,
                            flag_force_periodic);

         // Compute the local minima of each array we will exchange later
         if(x_i<x_min) x_min=x_i;
         if(y_i<y_min) y_min=y_i;
         if(z_i<z_min) z_min=z_i;
         if(h_i<h_min) h_min=h_i;
         if(f_i<f_min) f_min=f_i;
         if(v_i<f_min) v_min=v_i;
         if(w_i<f_min) w_min=w_i;

         // Populate z-array
         if(z_i>0){
            z_decomp[n_visible_local]=z_i;
            n_visible_local++;
         }
        
      }
      j_particle+=n_particles_species;
    }
  }

  // Compute sort indices  
  size_t *z_decomp_index;
  SID_Allreduce(&n_visible,&n_visible_local,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
  sort(z_decomp,n_visible_local,&z_decomp_index,SID_FLOAT,FALSE,SORT_COMPUTE_RANK,FALSE);
  
  // Set domain decomposition
  size_t i_particle_start;
  size_t i_particle_stop;
  size_t n_particles_local;
  size_t n_particles_alloc;
  i_particle=0;
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
     if(SID.My_rank==i_rank){
        i_particle_start=i_particle;
        if(SID.I_am_last_rank)
           i_particle_stop=n_visible-1;
        else
           i_particle_stop=(int)((float)(n_visible-1-i_particle_start)/(float)(SID.n_proc-i_rank));
        i_particle=i_particle_stop;
     }
     SID_Bcast(&i_particle,sizeof(size_t),i_rank,SID.COMM_WORLD);
  }
  n_particles_local=i_particle_stop-i_particle_start+1;
  SID_Allreduce(&n_particles_alloc,&n_particles_local,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);

  // Perform domain decomposition
  size_t  j_particle_start;
  size_t  j_particle_stop;
  size_t  n_particles_rank;
  float  *x_buffer;
  float  *y_buffer;
  float  *z_buffer;
  float  *h_buffer;
  float  *f_buffer;
  float  *v_buffer;
  float  *w_buffer;
  (*x)        =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  (*y)        =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  (*z)        =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  (*h_smooth) =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  (*f_stretch)=(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  x_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  y_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  z_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  h_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  f_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  v_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  w_buffer    =(float *)SID_malloc(sizeof(float)*(n_particles_rank));
  for(i_rank=0,k_particle=0;i_rank<SID.n_proc;i_rank++){
     // Communicate the particle ranges for the i_rank'th rank
     j_particle_start=i_particle_start;
     j_particle_stop =i_particle_stop;
     n_particles_rank=n_particles_local;
     SID_Bcast(&j_particle_start,sizeof(size_t),i_rank,SID.COMM_WORLD);
     SID_Bcast(&j_particle_stop, sizeof(size_t),i_rank,SID.COMM_WORLD);
     SID_Bcast(&n_particles_rank,sizeof(size_t),i_rank,SID.COMM_WORLD);
     // Initialize the exchange buffers
     for(i_particle=0;i_particle<n_particles_rank;i_particle++){
        x_buffer[i_particle]=x_min;
        y_buffer[i_particle]=y_min;
        z_buffer[i_particle]=z_min;
        h_buffer[i_particle]=0.;
        f_buffer[i_particle]=0.;
     }
     // Search local particles for those which belong on the i_rank'th rank
     for(i_particle=0;i_particle<n_particles_local;i_particle++){
        k_particle=z_decomp_index[i_particle]-j_particle_start;
        if(k_particle>=0 && k_particle<n_particles_rank){

           // Set the preoperties of the particle to be mapped
           set_particle_map_quantities(&mq,TRUE,i_particle,&x_i,&y_i,&z_i,&h_i,&v_i,&w_i);
   
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
                              flag_comoving,
                              flag_force_periodic);

           // Generate tretch factors
           switch(flag_plane_parallel){
              case FALSE:
                 if(z_i>0.)
                    f_i=d_o/z_i;
                 else
                    f_i=0.;
                 break;
              case TRUE:
                 f_i=1.;
                 break;
           }

           // Populate buffers
           x_buffer[k_particle]=x_i;
           y_buffer[k_particle]=y_i;
           z_buffer[k_particle]=z_i;
           h_buffer[k_particle]=h_i;
           f_buffer[k_particle]=f_i;
           v_buffer[k_particle]=v_i;
           w_buffer[k_particle]=w_i;
        }
     }
     // Perform buffer exchange
     SID_Reduce(x_buffer,(*x),        n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
     SID_Reduce(y_buffer,(*y),        n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
     SID_Reduce(z_buffer,(*z),        n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
     SID_Reduce(h_buffer,(*h_smooth), n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
     SID_Reduce(f_buffer,(*f_stretch),n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
     SID_Reduce(v_buffer,(*value),    n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
     SID_Reduce(w_buffer,(*weight),   n_particles_rank,SID_FLOAT,SID_MAX,i_rank,SID.COMM_WORLD);
  }
  SID_free(SID_FARG z_decomp_index);
  SID_free(SID_FARG x_buffer);
  SID_free(SID_FARG y_buffer);
  SID_free(SID_FARG z_buffer);
  SID_free(SID_FARG h_buffer);
  SID_free(SID_FARG f_buffer);
  SID_free(SID_FARG v_buffer);
  SID_free(SID_FARG w_buffer);

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
  float      *weight;
  float      *value;
  int        kx_min;
  int        kx_max;
  int        ky_min;
  int        ky_max;
  int        kz_min;
  int        kz_max;
  double    *numerator;
  double    *denominator;
  double    *z_image;
  int        py;
  int        px;
  int        pixel_pos;
  double     part_pos_x;
  double     part_pos_y;
  double     part_pos_z;
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
  int          flag_force_periodic;;
  double       expansion_factor;
  ADaPS       *transfer;
  double       h_Hubble;
  double       near_field;
  double       stereo_offset;
  int          i_x,i_y;
  int         *mask_buffer;

  int          i_image;
  int          camera_mode;

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

  if(check_mode_for_flag(camera_mode,CAMERA_STEREO))
    n_images=6;
  else
    n_images=2;

  for(i_image=0;i_image<n_images;i_image++){

    switch(i_image){
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
      case 2:
        parameter     =render->camera->RGB_param;
        transfer      =render->camera->RGB_transfer;
        image         =render->camera->image_RGB_left->values;
        z_image       =NULL;
        mask          =render->camera->mask_RGB_left;
        stereo_offset =-d_o/render->camera->stereo_ratio;
        FOV_x        +=2.*d_o/render->camera->stereo_ratio;
        break;
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
      case 4:
        parameter     =render->camera->RGB_param;
        transfer      =render->camera->RGB_transfer;
        image         =render->camera->image_RGB_right->values;
        z_image       =NULL;
        mask          =render->camera->mask_RGB_right;
        stereo_offset =d_o/render->camera->stereo_ratio;
        break;
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

    // Avoid overflows
    xmin  = -FOV_x/2.; // Things will be centred on (x_o,y_o,z_o) later
    ymin  = -FOV_y/2.; // Things will be centred on (x_o,y_o,z_o) later
  
    // Compute image scales
    n_pixels    =nx*ny;
    pixel_size_x=FOV_x/(double)nx;
    pixel_size_y=FOV_y/(double)ny;
    pixel_size  =0.5*(pixel_size_x+pixel_size_y);
    pixel_area  =pixel_size_x*pixel_size_y;
    if(fabs((pixel_size_x-pixel_size_y)/pixel_size_x)>1e-4)
      SID_log_warning("pixels are not square by %7.3f%%",0,fabs((pixel_size_x-pixel_size_y)/pixel_size_x)*1e2);

    // Initialize make_map
    init_make_map(plist,
                  parameter,
                  transfer,
                  x_o,y_o,z_o,
                  x_c,y_c,z_c,
                  box_size,FOV_x,
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
                  &n_particles);

    // Allocate and initialize image arrays
    if(flag_weigh)
      denominator=(double *)SID_malloc(sizeof(double)*n_pixels);
    else
      denominator=NULL;
    numerator=image;
    for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
      image[i_pixel]    =0.;
      numerator[i_pixel]=0.;
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

    // Generate the smoothing kernal
    set_sph_kernel(plist,SPH_KERNEL_GADGET|SPH_KERNEL_2D);
    kernel_radius      =(double *)ADaPS_fetch(plist->data, "sph_kernel_radius");
    kernel_table       =(double *)ADaPS_fetch(plist->data, "sph_kernel_2d");  
    kernel_table_avg   =((double *)ADaPS_fetch(plist->data,"sph_kernel_2d_dA"))[0];  
    radius_kernel_norm =kernel_radius[N_KERNEL_TABLE];
    radius_kernel_norm2=radius_kernel_norm*radius_kernel_norm;

    // Perform projection
    SID_log("Performing projection...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_particle=0,j_particle=0,k_particle++;i_particle<n_particles;i_particle++){
      part_h_z  =(double)h_smooth[i_particle];
      part_pos_z=(double)z[i_particle];
      if(part_pos_z>near_field){
        part_h_xy    =part_h_z*f_stretch[i_particle];
        part_pos_x   =(double)(x[i_particle]*f_stretch[i_particle]);
        part_pos_y   =(double)(y[i_particle]*f_stretch[i_particle]);
        px           =(int)((part_pos_x-xmin)/pixel_size_x+ONE_HALF);
        py           =(int)((part_pos_y-ymin)/pixel_size_y+ONE_HALF);
        pixel_pos    =py+px*ny;
        radius_kernel=radius_kernel_norm*part_h_xy;
        kx_min=(int)((part_pos_x-radius_kernel-xmin)/pixel_size_x);
        kx_max=(int)((part_pos_x+radius_kernel-xmin)/pixel_size_x+ONE_HALF);
        ky_min=(int)((part_pos_y-radius_kernel-ymin)/pixel_size_y);
        ky_max=(int)((part_pos_y+radius_kernel-ymin)/pixel_size_y+ONE_HALF);
        for(kx=kx_min,pixel_pos_x=xmin+(kx_min+0.5)*pixel_size_x;kx<=kx_max;kx++,pixel_pos_x+=pixel_size_x){
          if(kx>=0 && kx<nx){
            for(ky=ky_min,pixel_pos_y=ymin+(ky_min+0.5)*pixel_size_y;ky<=ky_max;ky++,pixel_pos_y+=pixel_size_y){
              if(ky>=0 && ky<ny){
                radius2_norm=
                  1./(part_h_xy*part_h_xy);
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
                    numerator[pos]  +=(double)value[i_particle]*(double)weight[i_particle]*kernel;
                    denominator[pos]+=(double)weight[i_particle]*kernel;
                    if(z_image!=NULL)
                      z_image[pos]+=(double)z[i_particle]*(double)weight[i_particle]*kernel;
                  }
                  else{
                    numerator[pos]+=(double)value[i_particle]*kernel;
                    if(z_image!=NULL)
                      z_image[pos]+=(double)z[i_particle]*(double)value[i_particle]*kernel;
                  }
                  mask[pos] =TRUE;
                }
              }
            }
          }
        }
      }
    }
    SID_Barrier(SID.COMM_WORLD);
    SID_log("Done.",SID_LOG_CLOSE);

    SID_log("Image normalization, etc...",SID_LOG_OPEN|SID_LOG_TIMER);

    // Add results from all ranks if this is being run in parallel
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

    // Normalize z-frame
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

    // Take log_10 if needed
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
    if(denominator!=NULL)
      SID_free(SID_FARG denominator);
    if(flag_weigh)
      SID_free(SID_FARG weight);
    SID_free(SID_FARG value);
  
    SID_log("Done.",SID_LOG_CLOSE);
  }

}
