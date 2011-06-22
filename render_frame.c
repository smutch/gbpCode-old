#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>
#include <gbpRender.h>

void rotate_particle(double  x_hat,
                     double  y_hat,
                     double  z_hat,
                     double  theta,
                     REAL   *x_i,
                     REAL   *y_i,
                     REAL   *z_i);
void rotate_particle(double  x_hat,
                     double  y_hat,
                     double  z_hat,
                     double  theta,
                     REAL   *x_i,
                     REAL   *y_i,
                     REAL   *z_i){
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

void init_make_map(plist_info  *plist,
                   char        *parameter,
                   ADaPS       *transfer_list,
                   double  x_o,
                   double  y_o,
                   double  z_o,
                   double  x_c,
                   double  y_c,
                   double  z_c,
                   double  box_size,
                   double  FOV,
                   int          mode,
                   double   expansion_factor,
                   double   near_field,
                   int          flag_comoving,
                   int         *flag_weigh,
                   int         *flag_line_integral,
                   REAL       **x,
                   REAL       **y,
                   REAL       **z,
                   REAL       **h_smooth,
                   REAL       **f_stretch,
                   REAL       **value,
                   REAL       **weight,
                   size_t      *n_particles);
void init_make_map(plist_info  *plist,
                   char        *parameter,
                   ADaPS       *transfer_list,
                   double   x_o,
                   double   y_o,
                   double   z_o,
                   double   x_c,
                   double   y_c,
                   double   z_c,
                   double   box_size,
                   double   FOV,
                   int      mode,
                   double   expansion_factor,
                   double   near_field,
                   int      flag_comoving,
                   int     *flag_weigh,
                   int     *flag_line_integral,
                   REAL   **x,
                   REAL   **y,
                   REAL   **z,
                   REAL   **h_smooth,
                   REAL   **f_stretch,
                   REAL   **value,
                   REAL   **weight,
                   size_t  *n_particles){
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
  REAL    *x_temp;
  REAL    *y_temp;
  REAL    *z_temp;
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
  
  SID_log("Initializing projection-space...",SID_LOG_OPEN|SID_LOG_TIMER);

  half_box_size=0.5*box_size;

  // Initialize array that will tell us what positions need to be used
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
    ptype_used[i_type]=FALSE;

  // Process parameter dependant stuff
  (*flag_weigh)        =FALSE;
  (*flag_line_integral)=FALSE;
  if(!strcmp(parameter,"Sigma_M_dark")){
    ptype_used[GADGET_TYPE_DARK]=TRUE;
    (*flag_weigh)               =FALSE;
    (*flag_line_integral)       =FALSE;
    (*n_particles)=((size_t *)ADaPS_fetch(plist->data,"n_dark"))[0];
    (*value)      = (REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
    (*weight)     = NULL;
    if(ADaPS_exist(plist->data,"rho_dark")){
      rho=(float *)ADaPS_fetch(plist->data,"rho_dark");
      for(i_particle=0;i_particle<(*n_particles);i_particle++)
        (*value)[i_particle]=(REAL)rho[i_particle];
    }
    else
      SID_trap_error("no denisities available to compute Sigma_M_dark in make_map",ERROR_LOGIC);
  }
  else if(!strcmp(parameter,"tau_dark")){
    ptype_used[GADGET_TYPE_DARK]=TRUE;
    (*flag_weigh)               =FALSE;
    (*flag_line_integral)       =TRUE;
    (*n_particles)=((size_t *)ADaPS_fetch(plist->data,"n_dark"))[0];
    (*value)      = (REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
    (*weight)     = NULL;
    if(ADaPS_exist(plist->data,"rho_dark")){
      rho=(float *)ADaPS_fetch(plist->data,"rho_dark");
      if(ADaPS_exist(transfer_list,"rho_dark")){
        if(ADaPS_exist(transfer_list,"rho_dark_log")) 
          flag_log=TRUE;
        transfer=(interp_info *)ADaPS_fetch(transfer_list,"rho_dark");
        for(i_particle=0;i_particle<(*n_particles);i_particle++){
          if(flag_log)
            transfer_val=MAX(0.,MIN(1.,interpolate(transfer,take_log10(rho[i_particle]))));
          else
            transfer_val=MAX(0.,MIN(1.,interpolate(transfer,rho[i_particle])));
          (*value)[i_particle]=(REAL)(rho[i_particle]*transfer_val);
        }
      }
      else{ 
        for(i_particle=0;i_particle<(*n_particles);i_particle++)
          (*value)[i_particle]=(REAL)rho[i_particle];
      }
      if(!flag_comoving){
        for(i_particle=0;i_particle<(*n_particles);i_particle++)
          ((*value)[i_particle])/=pow(expansion_factor,3.);
      }
    }
    else if(ADaPS_exist(plist->data,"mass_array_dark")){
      mass_array=((double *)ADaPS_fetch(plist->data,"mass_array_dark"))[0];
      for(i_particle=0;i_particle<(*n_particles);i_particle++)
        (*value)[i_particle]=(REAL)mass_array;
    }
    else
      SID_trap_error("no masses available to compute tau_dark in make_map",ERROR_LOGIC);
  }
  else if(!strcmp(parameter,"sigma_v_dark")){
    ptype_used[GADGET_TYPE_DARK]=TRUE;
    (*flag_weigh)               =TRUE;
    (*flag_line_integral)       =TRUE;
    (*n_particles)=((size_t *)ADaPS_fetch(plist->data,"n_dark"))[0];
    (*value)      = (REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
    (*weight)     = (REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
    // Compute a column-weighted average (using densities)
    if(ADaPS_exist(plist->data,"rho_dark")){
      rho=(float *)ADaPS_fetch(plist->data,"rho_dark");
      if(ADaPS_exist(transfer_list,"rho_dark")){
        if(ADaPS_exist(transfer_list,"rho_dark_log")) 
          flag_log=TRUE;
        transfer=(interp_info *)ADaPS_fetch(transfer_list,"rho_dark");
        for(i_particle=0;i_particle<(*n_particles);i_particle++){
          if(flag_log)
            transfer_val=MAX(0.,MIN(1.,interpolate(transfer,take_log10(rho[i_particle]))));
          else
            transfer_val=MAX(0.,MIN(1.,interpolate(transfer,rho[i_particle])));
          (*weight)[i_particle]=(REAL)(rho[i_particle]*transfer_val);
        }
      }
      else{
        for(i_particle=0;i_particle<(*n_particles);i_particle++)
          (*weight)[i_particle]=(REAL)rho[i_particle];
      }
      if(!flag_comoving){
        for(i_particle=0;i_particle<(*n_particles);i_particle++)
          ((*weight)[i_particle])/=pow(expansion_factor,3.);
      }
    }
    else if(ADaPS_exist(plist->data,"mass_array_dark")){
      mass_array=((double *)ADaPS_fetch(plist->data,"mass_array_dark"))[0];
      for(i_particle=0;i_particle<(*n_particles);i_particle++)
        (*weight)[i_particle]=(REAL)mass_array;
    }
    else
      SID_trap_error("no masses available to compute sigma_v_dark in make_map",ERROR_LOGIC);

    // Use sigma_v for values
    if(ADaPS_exist(plist->data,"sigma_v_dark")){
      sigma_v=(float *)ADaPS_fetch(plist->data,"sigma_v_dark");
      if(ADaPS_exist(transfer_list,"sigma_v_dark")){
        if(ADaPS_exist(transfer_list,"sigma_v_dark_log")) 
          flag_log=TRUE;
        transfer=(interp_info *)ADaPS_fetch(transfer_list,"sigma_v_dark");
        for(i_particle=0;i_particle<(*n_particles);i_particle++){
          if(flag_log)
            transfer_val=MAX(0.,MIN(1.,interpolate(transfer,take_log10(sigma_v[i_particle]))));
          else
            transfer_val=MAX(0.,MIN(1.,interpolate(transfer,sigma_v[i_particle])));
          (*value)[i_particle]=(REAL)(sigma_v[i_particle]*transfer_val);
        }
      } 
      else{
        for(i_particle=0;i_particle<(*n_particles);i_particle++)
          (*value)[i_particle]=(REAL)sigma_v[i_particle];
      }
    }
    else
      SID_trap_error("No sigma_v_dark's available for make_map",ERROR_LOGIC);
  }
  else
    SID_trap_error("Unknown quantity in sph_project {%s}",ERROR_LOGIC,parameter);

  // Set-up particle coordinates and smoothing lengths
  (*x)        =(REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
  (*y)        =(REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
  (*z)        =(REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
  (*h_smooth) =(REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
  (*f_stretch)=(REAL *)SID_malloc(sizeof(REAL)*(*n_particles));
  for(i_type=0,j_particle=0;i_type<N_GADGET_TYPE;i_type++){
    if(ptype_used[i_type] && ADaPS_exist(plist->data,"n_%s",plist->species[i_type])){
      n_particles_species=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_type]))[0];
      // Fetch coordinates
      x_temp=(REAL *)ADaPS_fetch(plist->data,"x_%s",plist->species[i_type]);
      y_temp=(REAL *)ADaPS_fetch(plist->data,"y_%s",plist->species[i_type]);
      z_temp=(REAL *)ADaPS_fetch(plist->data,"z_%s",plist->species[i_type]);
      for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++){
        (*x)[k_particle]=(REAL)(x_temp[i_particle]);
        (*y)[k_particle]=(REAL)(y_temp[i_particle]);
        (*z)[k_particle]=(REAL)(z_temp[i_particle]);
      }
      // Fetch smoothing lengths ...
      if(ADaPS_exist(plist->data,"r_smooth_%s",plist->species[i_type])){
        h_smooth_temp=(float   *)ADaPS_fetch(plist->data,"r_smooth_%s",plist->species[i_type]);
        for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++)
          (*h_smooth)[k_particle]=(REAL)(h_smooth_temp[i_particle]);
      }
      // ... if not present, create from densities ...
      else if(ADaPS_exist(plist->data,"rho_%s",plist->species[i_type])){
        rho=(float *)ADaPS_fetch(plist->data,"rho_%s",plist->species[i_type]);
        if(ADaPS_exist(plist->data,"M_dark")){
          mass=(double *)ADaPS_fetch(plist->data,"M_dark");
          for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++)
            (*h_smooth)[k_particle]=(REAL)pow(mass[i_particle]/(double)(rho[i_particle]),ONE_THIRD);        
        }
        else{
          if(ADaPS_exist(plist->data,"mass_array_dark"))
            mass_array=((double *)ADaPS_fetch(plist->data,"mass_array_dark"))[0];
          else
            mass_array=1.;
          for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++)
            (*h_smooth)[k_particle]=(REAL)pow(mass_array/(double)(rho[i_particle]),ONE_THIRD);        
        }
      }
      // ... if no densities are present, set to zero.
      else if(ADaPS_exist(plist->data,"r_smooth_%s",plist->species[i_type])){
        for(i_particle=0,k_particle=j_particle;i_particle<n_particles_species;i_particle++,k_particle++)
          (*h_smooth)[k_particle]=0.;        
      }
      j_particle+=n_particles_species;
    }
  }

  // Compute the angle and the axis of the rotation
  //   needed to place the camera at (0,0,-d_o)
  d_x_o=x_o-x_c;
  d_y_o=y_o-y_c;
  d_z_o=z_o-z_c;
  d_o  =sqrt(pow(d_x_o,2.)+
             pow(d_y_o,2.));
  x_hat     = d_y_o/d_o;
  y_hat     =-d_x_o/d_o;
  z_hat     = 0.;
  d_o  =sqrt(pow(d_x_o,2.)+
             pow(d_y_o,2.)+
             pow(d_z_o,2.));
  theta     =acos(d_z_o/d_o);
  if(sqrt(d_x_o*d_x_o+d_y_o*d_y_o)>0.){
    theta_roll=acos(-d_y_o/sqrt(d_x_o*d_x_o+d_y_o*d_y_o));
    if(d_x_o<0.)
      theta_roll=TWO_PI-theta_roll;
  }
  else
    theta_roll=0.;

  // Compute projection-space coordinates
  for(i_particle=0;i_particle<(*n_particles);i_particle++){
    
    // Shift origen to focus position
    (*x)[i_particle]-=x_o;
    (*y)[i_particle]-=y_o;
    (*z)[i_particle]-=z_o;

    if(!flag_comoving){
      force_periodic(&((*x)[i_particle]),-0.5*(REAL)box_size,(REAL)box_size);
      force_periodic(&((*y)[i_particle]),-0.5*(REAL)box_size,(REAL)box_size);
      force_periodic(&((*z)[i_particle]),-0.5*(REAL)box_size,(REAL)box_size);
    }

    // Rotate about origen to place camera position on z-axis
    rotate_particle((double)x_hat,
                    (double)y_hat,
                    (double)z_hat,
                    (double)theta,
                    &((*x)[i_particle]),
                    &((*y)[i_particle]),
                    &((*z)[i_particle]));

    // Apply roll angle
    rotate_particle((double)0.,
                    (double)0.,
                    (double)1.,
                    (double)(-theta_roll),
                    &((*x)[i_particle]),
                    &((*y)[i_particle]),
                    &((*z)[i_particle]));
                    
    // Convert to proper coordinates?
    if(!flag_comoving){
      ((*x)[i_particle])       *=expansion_factor;
      ((*y)[i_particle])       *=expansion_factor;
      ((*z)[i_particle])       *=expansion_factor;
      ((*h_smooth)[i_particle])*=expansion_factor;
    }

    // Shift image plane
    (*z)[i_particle]+=d_o;

    // Compute projection stretch
    if((*z)[i_particle]>0.)
      (*f_stretch)[i_particle]=d_o/(*z)[i_particle];
    else
      (*f_stretch)[i_particle]=0.;
  }

  // Apply a (optional) near-field correction
  if(near_field>0.){
    if((*flag_weigh)){
      for(i_particle=0;i_particle<(*n_particles);i_particle++){
        z_test=(*z)[i_particle];
        if(z_test>0. && z_test<(2.*near_field)){
          z_test/=near_field;
          (*weight)[i_particle]*=(1.-exp(-z_test*z_test));
        }
      }
    }
    else{
      for(i_particle=0;i_particle<(*n_particles);i_particle++){
        z_test=(*z)[i_particle];
        if(z_test>0. && z_test<(2.*near_field)){
          z_test/=near_field;
          (*value)[i_particle]*=(1.-exp(-z_test*z_test));
        }
      }
    }
  }

  SID_log("Done.",SID_LOG_CLOSE);
}

void render_frame(render_info  *render){
  size_t     i_particle;
  int        i_pixel;              
  size_t     i_kernel;
  size_t     n_pixels;
  size_t     n_unmasked;
  double     xmin,ymin;
  double     pixel_size;
  double     pixel_size_x;
  double     pixel_size_y;
  double     pixel_area;
  REAL    *weight;
  REAL    *value;
  int        kx_min;
  int        kx_max;
  int        ky_min;
  int        ky_max;
  int        kz_min;
  int        kz_max;
  int       *mask_local;
  double    *numerator;
  double    *numerator_local;
  double    *denominator;
  double    *denominator_local;
  double    *z_image;
  double    *z_frame;
  double    *z_frame_local;
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
  REAL      *x;
  REAL      *y;
  REAL      *z;
  REAL      *h_smooth;
  REAL      *f_stretch;
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

  double     z_avoidance;

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
  int         *mask;
  int          mode;
  char        *parameter;
  plist_info  *plist;
  int          flag_comoving;
  double       expansion_factor;
  ADaPS       *transfer;
  double       h_Hubble;
  double       near_field;

  int          i_image;

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
  flag_comoving   =render->flag_comoving;
  expansion_factor=render->camera->perspective->time;
  box_size        =((double *)ADaPS_fetch(plist->data,"box_size"))[0];
  h_Hubble        =render->h_Hubble;
  near_field      =render->near_field*d_o;
  if(nx>=ny){
    FOV_y=render->camera->perspective->FOV;
    FOV_x=FOV_y*(double)nx/(double)ny;    
  }
  else{
    FOV_x=render->camera->perspective->FOV;
    FOV_y=FOV_x*(double)nx/(double)ny;        
  }

  for(i_image=0;i_image<2;i_image++){

    switch(i_image){
      case 0:
        parameter=render->camera->RGB_param;
        image    =render->camera->image_RGB->values;
        mask     =render->camera->mask_RGB;
        transfer =render->camera->RGB_transfer;
        z_image  =NULL;
        break;
      case 1:
        parameter=render->camera->Y_param;
        image    =render->camera->image_Y->values;
        mask     =render->camera->mask_Y;
        transfer =render->camera->Y_transfer;
        z_image  =render->camera->image_Z->values;
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
                  flag_comoving,
                  &flag_weigh,
                  &flag_line_integral,
                  &x,&y,&z,
                  &h_smooth,
                  &f_stretch,
                  &value,
                  &weight,
                  &n_particles);

    // Allocate and initialize image arrays
    numerator_local=(double *)SID_malloc(sizeof(double)*n_pixels);
    if(flag_weigh)
      denominator_local=(double *)SID_malloc(sizeof(double)*n_pixels);
    if(z_frame!=NULL)
      z_frame_local=(double *)SID_malloc(sizeof(double)*n_pixels);
    mask_local=(int   *)SID_malloc(sizeof(int)*n_pixels);
    for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
      image[i_pixel]          =0.;
      numerator_local[i_pixel]=0.;
      if(flag_weigh)
        denominator_local[i_pixel]=0.;
      if(z_image!=NULL)
        z_frame_local[i_pixel] =0.;
      mask_local[i_pixel]=FALSE;
    }

    d_x_o=x_o-x_c;
    d_y_o=y_o-y_c;
    d_z_o=z_o-z_c;
    d_o  =sqrt(pow(d_x_o,2.)+
               pow(d_y_o,2.)+
               pow(d_z_o,2.));

    z_avoidance=0.05*d_o;

    // Generate the smoothing kernal
    set_sph_kernel(plist,SPH_KERNEL_GADGET|SPH_KERNEL_2D);
    kernel_radius      =(double *)ADaPS_fetch(plist->data, "sph_kernel_radius");
    kernel_table       =(double *)ADaPS_fetch(plist->data, "sph_kernel_2d");  
    kernel_table_avg   =((double *)ADaPS_fetch(plist->data,"sph_kernel_2d_dA"))[0];  
    radius_kernel_norm =kernel_radius[N_KERNEL_TABLE];
    radius_kernel_norm2=radius_kernel_norm*radius_kernel_norm;

    // Perform projection
    SID_log("Performing projection...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(i_particle=0;i_particle<n_particles;i_particle++){
      part_h_z  =(double)h_smooth[i_particle];
      part_pos_z=(double)z[i_particle];
      if(part_pos_z>z_avoidance){
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
                  if(flag_weigh){
                    numerator_local[pos]  +=(double)value[i_particle]*(double)weight[i_particle]*kernel;
                    denominator_local[pos]+=(double)weight[i_particle]*kernel;
                    if(z_image!=NULL)
                      z_frame_local[pos]+=(double)z[i_particle]*(double)weight[i_particle]*kernel;
                  }
                  else{
                    numerator_local[pos]+=(double)value[i_particle]*kernel;
                    if(z_image!=NULL)
                      z_frame_local[pos]+=(double)z[i_particle]*(double)value[i_particle]*kernel;
                  }
                  mask_local[pos] =TRUE;
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
#ifdef USE_MPI
    MPI_Allreduce(mask_local,mask,n_pixels,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    SID_free((void **)&mask_local);
    numerator=(double *)SID_malloc(sizeof(double)*n_pixels);
    MPI_Allreduce(numerator_local,numerator,n_pixels,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(flag_weigh){
      denominator=(double *)SID_malloc(sizeof(double)*n_pixels);
      MPI_Allreduce(denominator_local,denominator,n_pixels,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    if(z_image!=NULL){
      z_frame=(double *)SID_malloc(sizeof(double)*n_pixels);
      MPI_Allreduce(z_frame_local,z_frame,n_pixels,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
#else
    mask       =mask_local;
    numerator  =numerator_local;
    denominator=denominator_local;
    z_frame    =z_frame_local;
#endif

    // Normalize image (if needed)
    if(flag_weigh){
      for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
        if(mask[i_pixel])
          image[i_pixel]=numerator[i_pixel]/denominator[i_pixel];
      }
    }
    else{
      for(i_pixel=0;i_pixel<n_pixels;i_pixel++)
        image[i_pixel]=numerator[i_pixel];
    }

    // Normalize z-frame
    if(z_image!=NULL){
      if(flag_weigh){
        for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
            z_image[i_pixel]=z_frame[i_pixel]/denominator[i_pixel];
        }
      }
      else{
        for(i_pixel=0;i_pixel<n_pixels;i_pixel++){
          if(mask[i_pixel])
            z_image[i_pixel]=z_frame[i_pixel]/numerator[i_pixel];
        }
      }
    }

#ifdef USE_MPI
    SID_free(SID_FARG numerator);
    if(flag_weigh)
      SID_free(SID_FARG denominator);
    if(z_image!=NULL)
      SID_free(SID_FARG z_frame);
#endif

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
    for(i_pixel=0,n_unmasked=0;i_pixel<n_pixels;i_pixel++) if(mask[i_pixel]) n_unmasked++;
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
    SID_free((void **)&x);
    SID_free((void **)&y);
    SID_free((void **)&z);
    SID_free((void **)&h_smooth);
    SID_free((void **)&f_stretch);
    SID_free((void **)&value);
    SID_free((void **)&numerator_local);
    if(flag_weigh){
      SID_free((void **)&weight);
      SID_free((void **)&denominator_local);
    }
    if(z_image!=NULL)
      SID_free(SID_FARG z_frame_local);
  
    SID_log("Done.",SID_LOG_CLOSE);
  }
}
