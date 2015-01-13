#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>

#define FLAG_INTERP_MIN_BIN 4
#define FLAG_INTERP_MAX_BIN 8

int compute_group_analysis(halo_properties_info  *properties,
                           halo_profile_info     *profile,
                           double (*p_i_fctn) (void *,int,int), 
                           double (*v_i_fctn) (void *,int,int), 
                           size_t (*id_i_fctn)(void *,int), 
                           void   *params,
                           //size_t                *id_array,
                           //GBPREAL               *x_array,
                           //GBPREAL               *y_array,
                           //GBPREAL               *z_array,
                           //GBPREAL               *vx_array,
                           //GBPREAL               *vy_array,
                           //GBPREAL               *vz_array,
                           //size_t                *ids_sort_index,
                           double                 box_size,
                           double                 particle_mass,
                           int                    n_particles,
                           double                 expansion_factor,
                           int                    i_halo,
                           double                *x,
                           double                *y,
                           double                *z,
                           double                *vx,
                           double                *vy,
                           double                *vz,
                           double                *R,
                           size_t               **R_index_in,
                           int                    flag_manual_centre,
                           cosmo_info            *cosmo){
  size_t      *R_index;
  int          i,j;
  int          i_profile;
  int          j_profile;
  int          k_profile;
  int          n_profile;
  size_t       i_particle;
  size_t       j_particle;
  size_t       k_particle;
  int          next_bin_particle;
  interp_info *V_R_interpolate;
  interp_info *vir_interpolate;
  int          i_bin;
  int          n_in_bin;
  double       n_per_bin;
  int          n_cumulative;
  double       V1;
  double       V2;
  double       dV;
  double       dM;
  double       sigma_r_mean;
  double       sigma_t_mean;
  double       sigma_T_mean;
  double       sigma_P_mean;
  double       sigma_mean;
  double       x_COM_accumulator;
  double       y_COM_accumulator;
  double       z_COM_accumulator;
  double       vx_COM_accumulator;
  double       vy_COM_accumulator;
  double       vz_COM_accumulator;
  double       spin_x_accumulator;
  double       spin_y_accumulator;
  double       spin_z_accumulator;
  double       r_xy;
  double       v_tot,v_rad,v_tan;
  double       v_x_mean,v_y_mean,v_z_mean,v_rad_mean;
  double       shape_eigen_values[3];
  double       shape_eigen_vectors[3][3];
  double       x_COM,y_COM,z_COM,R_COM;
  double       r_c[MAX_PROFILE_BINS_P1];
  double       v_c[MAX_PROFILE_BINS_P1];
  double       r_interp[MAX_PROFILE_BINS];
  double       y_interp[MAX_PROFILE_BINS];
  size_t       n_bins_temp;
  double       V_max,R_max;
  double       Delta,Omega;
  double       norm;
  int          flag_interpolated=FALSE;
  const gsl_interp_type  *interp_type;
  double sigma_cor,sigma_halo;
  double M_cor,M_halo;
  double x_vir,gamma;

  double h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  double Omega_M =((double *)ADaPS_fetch(cosmo,"Omega_M"))[0];
  double redshift=z_of_a(expansion_factor);

  Delta=Delta_vir(redshift,cosmo);
  Omega=1.;

  // Initialize properties
  properties->id_MBP         =id_i_fctn(params,0);//id_array[index_MBP];
  properties->n_particles    =n_particles;
  properties->position_COM[0]=0.;
  properties->position_COM[1]=0.;
  properties->position_COM[2]=0.;
  properties->position_MBP[0]=(float)p_i_fctn(params,0,0);//(x_array[index_MBP]);
  properties->position_MBP[1]=(float)p_i_fctn(params,1,0);//(y_array[index_MBP]);
  properties->position_MBP[2]=(float)p_i_fctn(params,2,0);//(z_array[index_MBP]);
  properties->velocity_COM[0]=0.;
  properties->velocity_COM[1]=0.;
  properties->velocity_COM[2]=0.;
  properties->velocity_MBP[0]=(float)v_i_fctn(params,0,0);//(vx_array[index_MBP]);
  properties->velocity_MBP[1]=(float)v_i_fctn(params,1,0);//(vy_array[index_MBP]);
  properties->velocity_MBP[2]=(float)v_i_fctn(params,2,0);//(vz_array[index_MBP]);
  properties->M_vir          =0.;
  properties->R_vir          =0.;
  properties->R_halo         =0.;
  properties->R_max          =0.;
  properties->V_max          =0.;
  properties->sigma_v        =0.;
  properties->spin[0]        =0.;
  properties->spin[1]        =0.;
  properties->spin[2]        =0.;
  properties->q_triaxial     =0.;
  properties->s_triaxial     =0.;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      properties->shape_eigen_vectors[i][j]=0.;
  }

  // Set the number of profile bins and the number of particles per bin
  profile->n_bins=MAX(MIN_PROFILE_BINS,MIN((int)((6.2*log10((double)n_particles)-3.5)+((double)n_particles/1000.)+1),MAX_PROFILE_BINS));
  n_per_bin      =(double)(n_particles)/(double)profile->n_bins;

  // There's nothing to do if there are no particles
  if(n_particles==0)
    profile->n_bins=0;
  else{
    // Create a v_c(0)=0 bin
    r_c[0]=0.;
    v_c[0]=0.;

    // Initialize profiles
    for(i_bin=0;i_bin<profile->n_bins;i_bin++){
      profile->bins[i_bin].r_med           =0.;
      profile->bins[i_bin].r_max           =0.;
      profile->bins[i_bin].n_particles     =0;
      profile->bins[i_bin].M_r             =0.;
      profile->bins[i_bin].rho             =0.;
      profile->bins[i_bin].overdensity     =0.;
      profile->bins[i_bin].position_COM[0] =0.;
      profile->bins[i_bin].position_COM[1] =0.;
      profile->bins[i_bin].position_COM[2] =0.;
      profile->bins[i_bin].velocity_COM[0] =0.;
      profile->bins[i_bin].velocity_COM[1] =0.;
      profile->bins[i_bin].velocity_COM[2] =0.;
      profile->bins[i_bin].sigma_rad       =0.;
      profile->bins[i_bin].sigma_tan       =0.;
      profile->bins[i_bin].sigma_tot       =0.;
      profile->bins[i_bin].spin[0]         =0.;
      profile->bins[i_bin].spin[1]         =0.;
      profile->bins[i_bin].spin[2]         =0.;
      profile->bins[i_bin].q_triaxial      =0.;
      profile->bins[i_bin].s_triaxial      =0.;
      for(i=0;i<3;i++){
        for(j=0;j<3;j++)
          profile->bins[i_bin].shape_eigen_vectors[i][j]=0.;
      }
    }

    // Fill temporary arrays for particle positions, radii (all w.r.t MBP) and velocities
    //   Also, enforce periodic box on particle positions
    double x_cen;
    double y_cen;
    double z_cen;
    x_cen=(double)properties->position_MBP[0];
    y_cen=(double)properties->position_MBP[1];
    z_cen=(double)properties->position_MBP[2];
    if(flag_manual_centre){
       double x_cen_manual;
       double y_cen_manual;
       double z_cen_manual;
       // Compute a rough comoving centre
       for(j_particle=0;j_particle<n_particles;j_particle++){
         x[j_particle]=d_periodic(p_i_fctn(params,0,j_particle)-x_cen,box_size);//(double)(x_array[k_particle])-x_cen,box_size);
         y[j_particle]=d_periodic(p_i_fctn(params,1,j_particle)-y_cen,box_size);//(double)(y_array[k_particle])-y_cen,box_size);
         z[j_particle]=d_periodic(p_i_fctn(params,2,j_particle)-z_cen,box_size);//(double)(z_array[k_particle])-z_cen,box_size);
       }
       // Refine it with shrinking spheres
       int n_iterations;
       n_iterations=compute_centroid3D(NULL,
                                       x,
                                       y,
                                       z,
                                       n_particles,
                                       1e-3, // 1 kpc
                                       0.75,
                                       30,
                                       CENTROID3D_MODE_FACTOR|CENTROID3D_MODE_INPLACE,
                                       &x_cen_manual,
                                       &y_cen_manual,
                                       &z_cen_manual);
       x_cen+=x_cen_manual;
       y_cen+=y_cen_manual;
       z_cen+=z_cen_manual;
       properties->position_MBP[0]=x_cen;
       properties->position_MBP[1]=y_cen;
       properties->position_MBP[2]=z_cen;
       if(properties->position_MBP[0]< box_size) properties->position_MBP[0]+=box_size;
       if(properties->position_MBP[1]< box_size) properties->position_MBP[1]+=box_size;
       if(properties->position_MBP[2]< box_size) properties->position_MBP[2]+=box_size;
       if(properties->position_MBP[0]>=box_size) properties->position_MBP[0]-=box_size;
       if(properties->position_MBP[1]>=box_size) properties->position_MBP[1]-=box_size;
       if(properties->position_MBP[2]>=box_size) properties->position_MBP[2]-=box_size;
    }
    
    for(j_particle=0;j_particle<n_particles;j_particle++){
      // ... halo-centric particle positions ...
      x[j_particle]=expansion_factor*d_periodic(p_i_fctn(params,0,j_particle)-x_cen,box_size);//((double)x_array[k_particle])-x_cen,box_size);
      y[j_particle]=expansion_factor*d_periodic(p_i_fctn(params,1,j_particle)-y_cen,box_size);//((double)y_array[k_particle])-y_cen,box_size);
      z[j_particle]=expansion_factor*d_periodic(p_i_fctn(params,2,j_particle)-z_cen,box_size);//((double)z_array[k_particle])-z_cen,box_size);

      // ... velocities ...
      vx[j_particle]=v_i_fctn(params,0,j_particle);//(double)(vx_array[k_particle]);
      vy[j_particle]=v_i_fctn(params,1,j_particle);//(double)(vy_array[k_particle]);
      vz[j_particle]=v_i_fctn(params,2,j_particle);//(double)(vz_array[k_particle]);

      // ... particle radii ...
      R[j_particle]=sqrt(x[j_particle]*x[j_particle]+y[j_particle]*y[j_particle]+z[j_particle]*z[j_particle]);
    }
    
    // Sort particles by radius
    merge_sort((void *)R,(size_t)n_particles,R_index_in,SID_DOUBLE,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    R_index=(*R_index_in);

    // Use the average of the central 30 particles for the MBP velocity if we are
    //    manually computing centres
    if(flag_manual_centre){
       double vx_cen_temp=0.;
       double vy_cen_temp=0.;
       double vz_cen_temp=0.;
       int    n_cen =0;
       for(i_particle=0;i_particle<MIN(30,n_particles);i_particle++){
          vx_cen_temp+=vx[i_particle];
          vy_cen_temp+=vy[i_particle];
          vz_cen_temp+=vz[i_particle];
          n_cen++;
       }
       properties->velocity_MBP[0]=vx_cen_temp/(double)n_cen;
       properties->velocity_MBP[1]=vy_cen_temp/(double)n_cen;
       properties->velocity_MBP[2]=vz_cen_temp/(double)n_cen;
    }

    // We need the COM velocity at R_vir before we can get halo centric velocities.  Thus,
    //   we need the overdensity profile first
    x_COM_accumulator  =0.;
    y_COM_accumulator  =0.;
    z_COM_accumulator  =0.;
    vx_COM_accumulator =0.;
    vy_COM_accumulator =0.;
    vz_COM_accumulator =0.;
    V2                 =0.;
    n_cumulative       =0;
    for(i_bin=0,i_particle=0;i_bin<profile->n_bins;i_bin++,i_particle+=n_in_bin){

      V1=V2; // Volumes

      // ... particle numbers ...
      if(i_bin<profile->n_bins-1)
        n_in_bin=(int)((double)(i_bin+1)*n_per_bin)-i_particle;
      else
        n_in_bin=n_particles-i_particle;
      n_cumulative                    +=n_in_bin;
      profile->bins[i_bin].n_particles =n_in_bin;

      // ... mass profile ...
      profile->bins[i_bin].M_r=particle_mass*(double)n_cumulative;

      // ... binning radii ...
      if(n_in_bin%2==1)
        profile->bins[i_bin].r_med=(float)R[R_index[i_particle+n_in_bin/2]];
      else
        profile->bins[i_bin].r_med=0.5*(float)(R[R_index[i_particle+n_in_bin/2-1]]+R[R_index[i_particle+n_in_bin/2]]);
      profile->bins[i_bin].r_max=(float)R[R_index[i_particle+n_in_bin-1]];

      // ... COM positions and velocities ...
      for(j_particle=0;j_particle<n_in_bin;j_particle++){
        k_particle=R_index[i_particle+j_particle];
        x_COM_accumulator += x[k_particle];
        y_COM_accumulator += y[k_particle];
        z_COM_accumulator += z[k_particle];
        vx_COM_accumulator+=vx[k_particle];
        vy_COM_accumulator+=vy[k_particle];
        vz_COM_accumulator+=vz[k_particle];
      }
      profile->bins[i_bin].position_COM[0]=(float)(x_COM_accumulator/(double)n_cumulative);
      profile->bins[i_bin].position_COM[1]=(float)(y_COM_accumulator/(double)n_cumulative);
      profile->bins[i_bin].position_COM[2]=(float)(z_COM_accumulator/(double)n_cumulative);
      profile->bins[i_bin].velocity_COM[0]=(float)(vx_COM_accumulator/(double)n_cumulative);
      profile->bins[i_bin].velocity_COM[1]=(float)(vy_COM_accumulator/(double)n_cumulative);
      profile->bins[i_bin].velocity_COM[2]=(float)(vz_COM_accumulator/(double)n_cumulative);

      // ... density ...
      V2=FOUR_THIRDS_PI*profile->bins[i_bin].r_max*profile->bins[i_bin].r_max*profile->bins[i_bin].r_max; // Volume
      dV=V2-V1;
      dM=particle_mass*(double)n_in_bin;
      profile->bins[i_bin].rho        =(float)(dM/dV);
      profile->bins[i_bin].overdensity=(float)(profile->bins[i_bin].M_r/(V2*Omega*rho_crit_z(redshift,cosmo)));

      /// ... triaxiality ...
      compute_triaxiality(x,
                          y,
                          z,
                          (double)profile->bins[i_bin].position_COM[0],
                          (double)profile->bins[i_bin].position_COM[1],
                          (double)profile->bins[i_bin].position_COM[2],
                          box_size,
                          n_cumulative,
                          R_index,
                          shape_eigen_values,
                          shape_eigen_vectors);
      profile->bins[i_bin].q_triaxial=(float)(shape_eigen_values[1]/shape_eigen_values[0]);
      profile->bins[i_bin].s_triaxial=(float)(shape_eigen_values[2]/shape_eigen_values[0]);
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)
          profile->bins[i_bin].shape_eigen_vectors[i][j]=(float)shape_eigen_vectors[i][j];
    }

    // Interpolate to get R_vir
    flag_interpolated=FALSE;
    properties->R_halo=profile->bins[profile->n_bins-1].r_max;
    if(profile->n_bins>1){

      // Remove any small-radius monotonic increases from the interpolation interval
      j_profile=0;
      while(profile->bins[j_profile].overdensity<=profile->bins[j_profile+1].overdensity && j_profile<profile->n_bins-2)
        j_profile++;

      // Only keep decreasing bins
      n_bins_temp=0;
      r_interp[n_bins_temp]=take_log10((double)profile->bins[j_profile].r_max);
      y_interp[n_bins_temp]=take_log10((double)profile->bins[j_profile].overdensity);
      n_bins_temp++;
      for(i_profile=j_profile+1;i_profile<profile->n_bins;i_profile++){
        if(take_log10((double)profile->bins[i_profile].overdensity)<y_interp[n_bins_temp-1]){
          r_interp[n_bins_temp]=take_log10((double)profile->bins[i_profile].r_max);
          y_interp[n_bins_temp]=take_log10((double)profile->bins[i_profile].overdensity);
          n_bins_temp++;
        }
      }

      if(n_bins_temp>1){
        // Perform interpolation
        if(y_interp[0]>=take_log10(Delta) && y_interp[n_bins_temp-1]<=take_log10(Delta)){
          if(n_bins_temp>9)
            interp_type=gsl_interp_cspline;
          else
            interp_type=gsl_interp_linear;
          interp_type=gsl_interp_linear;
          init_interpolate(y_interp,r_interp,n_bins_temp,interp_type,&vir_interpolate);
          properties->R_vir       =(float)take_alog10(interpolate(vir_interpolate,take_log10(Delta)));
          free_interpolate(SID_FARG vir_interpolate,NULL);
          flag_interpolated=TRUE;
        }
        else if(y_interp[0]<take_log10(Delta)){
          properties->R_vir=(float)take_alog10(r_interp[0]);
          flag_interpolated=FLAG_INTERP_MIN_BIN;
        }
        else{
          properties->R_vir=(float)take_alog10(r_interp[n_bins_temp-1]);
          flag_interpolated=FLAG_INTERP_MAX_BIN;
        }
      }
      else{
        properties->R_vir=(float)take_alog10(r_interp[0]);
        flag_interpolated=FLAG_INTERP_MIN_BIN;
      }
    }
    else{
      properties->R_vir=profile->bins[0].r_max; 
      flag_interpolated=FLAG_INTERP_MIN_BIN;
    }

    // Set the interpolation method
    if(n_bins_temp>9)
      interp_type=gsl_interp_cspline;
    else
      interp_type=gsl_interp_linear;
    interp_type=gsl_interp_linear;

    // Compute v_COM(R_vir)
    for(i_profile=0;i_profile<profile->n_bins;i_profile++)
      r_interp[i_profile]=(double)profile->bins[i_profile].r_max;
    if(flag_interpolated==TRUE){
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].velocity_COM[0];
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->velocity_COM[0]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].velocity_COM[1];
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->velocity_COM[1]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].velocity_COM[2];
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->velocity_COM[2]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
    }
    else{
      if(flag_interpolated==FLAG_INTERP_MIN_BIN)      i_profile=0;
      else if(flag_interpolated==FLAG_INTERP_MAX_BIN) i_profile=profile->n_bins-1;
      else SID_trap_error("Unrecognized value for flag_interpolated {%d}.",flag_interpolated);
      properties->velocity_COM[0]=(float)profile->bins[i_profile].velocity_COM[0];
      properties->velocity_COM[1]=(float)profile->bins[i_profile].velocity_COM[1];
      properties->velocity_COM[2]=(float)profile->bins[i_profile].velocity_COM[2];
    }

    // Compute halo-centric particle velocities
    //   Subtract COM mean and add Hubble flow
    for(j_particle=0;j_particle<n_particles;j_particle++){
      vx[j_particle]+=x[j_particle]*H_convert(H_z(redshift,cosmo))-(double)properties->velocity_COM[0];
      vy[j_particle]+=y[j_particle]*H_convert(H_z(redshift,cosmo))-(double)properties->velocity_COM[1];
      vz[j_particle]+=z[j_particle]*H_convert(H_z(redshift,cosmo))-(double)properties->velocity_COM[2];
    }

    // Compute remaining profiles ...
    spin_x_accumulator=0.;
    spin_y_accumulator=0.;
    spin_z_accumulator=0.;
    V2                =0.;
    n_cumulative      =0;
    for(i_bin=0,i_particle=0;i_bin<profile->n_bins;i_bin++,i_particle+=n_in_bin){

      V1=V2; // Volumes

      // ... particle numbers ...
      if(i_bin<profile->n_bins-1)
        n_in_bin=(int)((double)(i_bin+1)*n_per_bin)-i_particle;
      else
        n_in_bin=n_particles-i_particle;
      n_cumulative+=n_in_bin;

      // ... spins and mean velocities ...
      v_x_mean  =0.;
      v_y_mean  =0.;
      v_z_mean  =0.;
      v_rad_mean=0.;
      for(j_particle=0;j_particle<n_in_bin;j_particle++){
        k_particle=R_index[i_particle+j_particle];
      
        // ... spins ...
        spin_x_accumulator+=(double)(y[k_particle]*vz[k_particle]-z[k_particle]*vy[k_particle]);
        spin_y_accumulator+=(double)(z[k_particle]*vx[k_particle]-x[k_particle]*vz[k_particle]);
        spin_z_accumulator+=(double)(x[k_particle]*vy[k_particle]-y[k_particle]*vx[k_particle]);

        // ... mean velocities (needed below for velocity dispersions) ...
        if(R[k_particle]>0.){
          v_rad      =(x[k_particle]*vx[k_particle]+y[k_particle]*vy[k_particle]+z[k_particle]*vz[k_particle])/R[k_particle];
          v_x_mean  +=vx[k_particle];
          v_y_mean  +=vy[k_particle];
          v_z_mean  +=vz[k_particle];
          v_rad_mean+=v_rad;
        }
      }
      profile->bins[i_bin].spin[0]=(float)(spin_x_accumulator)/n_cumulative;
      profile->bins[i_bin].spin[1]=(float)(spin_y_accumulator)/n_cumulative;
      profile->bins[i_bin].spin[2]=(float)(spin_z_accumulator)/n_cumulative;
      v_x_mean  /=(double)n_in_bin;
      v_y_mean  /=(double)n_in_bin;
      v_z_mean  /=(double)n_in_bin;
      v_rad_mean/=(double)n_in_bin;

      // ... velocity dispersions ...
      for(j_particle=0;j_particle<n_in_bin;j_particle++){
        k_particle=R_index[i_particle+j_particle];
        if(R[k_particle]>0.){
          v_tot=sqrt(pow(vx[k_particle]-v_x_mean,2.)+pow(vy[k_particle]-v_y_mean,2.)+pow(vz[k_particle]-v_z_mean,2.));
          v_rad=(x[k_particle]*(vx[k_particle]-v_x_mean)+y[k_particle]*(vy[k_particle]-v_y_mean)+z[k_particle]*(vz[k_particle]-v_z_mean))/R[k_particle];
          v_tan=sqrt(v_tot*v_tot-v_rad*v_rad);
          profile->bins[i_bin].sigma_tot+=(float)((v_tot)*(v_tot));
          profile->bins[i_bin].sigma_rad+=(float)((v_rad-v_rad_mean)*(v_rad-v_rad_mean));
          profile->bins[i_bin].sigma_tan+=(float)((v_tan)*(v_tan));
        }
      }
      profile->bins[i_bin].sigma_tot=(float)sqrt((double)profile->bins[i_bin].sigma_tot/(double)n_in_bin);
      profile->bins[i_bin].sigma_rad=(float)sqrt((double)profile->bins[i_bin].sigma_rad/(double)n_in_bin);
      profile->bins[i_bin].sigma_tan=(float)sqrt((double)profile->bins[i_bin].sigma_tan/(double)n_in_bin);
    
      // ... circular velocity; v_c(R) ...
      r_c[i_bin+1]=profile->bins[i_bin].r_max;
      v_c[i_bin+1]=sqrt(G_NEWTON*profile->bins[i_bin].M_r/(double)profile->bins[i_bin].r_max);
    }
    
    // Determine R_max and V_max from v_c(R)...
    R_max=(double)r_c[1]; // default for a monotonically increasing V_c(r)
    V_max=(double)v_c[1]; // default for a monotonically increasing V_c(r)
    if(profile->n_bins>1){

      // Remove any large-radius monotonic increases from the interpolation interval
      k_profile=profile->n_bins;
      while(v_c[k_profile-1]<=v_c[k_profile] && k_profile>1)
        k_profile--;
      if(v_c[0]<=v_c[1] && k_profile==1)
        k_profile--;

      // If the profile is not monotonically increasing ...
      if(k_profile>0){
        n_bins_temp=k_profile+1;
      
        // ...find the maximum (call its index j_profile)
        for(i_profile=0,j_profile=0;i_profile<n_bins_temp;i_profile++){
          if(v_c[i_profile]>v_c[j_profile])
            j_profile=i_profile;
        }

        // ...find bottom of range in which to search for maximum (call its index i_profile)
        i_profile=j_profile-1;
        while(v_c[i_profile]>=v_c[j_profile] && i_profile>0)
          i_profile--;
    
        // ...find top of range in which to search for maximum (call its index k_profile)
        k_profile=j_profile+1;
        while(v_c[k_profile]>=v_c[j_profile] && k_profile<n_bins_temp-1)
          k_profile++;

        //   ... perform interpolation
        V_max=(double)v_c[j_profile];
        R_max=(double)r_c[j_profile];
        if(i_profile<j_profile && j_profile<k_profile){
          if(n_bins_temp>9)
            interp_type=gsl_interp_cspline;
          else
            interp_type=gsl_interp_linear;
          interp_type=gsl_interp_linear;
          init_interpolate(r_c,v_c,n_bins_temp,gsl_interp_cspline,&V_R_interpolate);
          interpolate_maximum(V_R_interpolate,
                              r_c[i_profile],
                              r_c[j_profile],
                              r_c[k_profile],
                              0.05,
                              &R_max,
                              &V_max);
          free_interpolate(SID_FARG V_R_interpolate,NULL);
        }
      }
    }
    properties->R_max=(float)R_max;
    properties->V_max=(float)V_max;

    if(profile->n_bins>9)
      interp_type=gsl_interp_cspline;
    else
      interp_type=gsl_interp_linear;
    interp_type=gsl_interp_linear;

    // Perform normal interpolation from profiles to get the rest of the global quantities
    if(flag_interpolated==TRUE){
      //  ... COM positions ...
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].position_COM[0]/expansion_factor;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->position_COM[0]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].position_COM[1]/expansion_factor;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->position_COM[1]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].position_COM[2]/expansion_factor;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->position_COM[2]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      //  ... M_vir ...
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].M_r;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->M_vir=interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      //  ... sigma_v ...
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].sigma_tot;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->sigma_v=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      //   ... spin ...
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].spin[0];
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->spin[0]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].spin[1];
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->spin[1]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].spin[2];
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->spin[2]=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      //  ... triaxial axes ratios ...
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].q_triaxial;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->q_triaxial=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      for(i_profile=0;i_profile<profile->n_bins;i_profile++)
        y_interp[i_profile]=(double)profile->bins[i_profile].s_triaxial;
      init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
      properties->s_triaxial=(float)interpolate(vir_interpolate,properties->R_vir);
      free_interpolate(SID_FARG vir_interpolate,NULL);
      // ... shape eigen vectors ...
      for(i=0;i<3;i++){
        for(j=0;j<3;j++){
          for(i_profile=0;i_profile<profile->n_bins;i_profile++)
            y_interp[i_profile]=(double)cos(profile->bins[i_profile].shape_eigen_vectors[i][j]);
          init_interpolate(r_interp,y_interp,profile->n_bins,interp_type,&vir_interpolate);
          properties->shape_eigen_vectors[i][j]=(float)acos(MAX(0,MIN(1.,interpolate(vir_interpolate,properties->R_vir))));
          free_interpolate(SID_FARG vir_interpolate,NULL);
        }
        norm=sqrt(properties->shape_eigen_vectors[i][0]*properties->shape_eigen_vectors[i][0]+
                  properties->shape_eigen_vectors[i][1]*properties->shape_eigen_vectors[i][1]+
                  properties->shape_eigen_vectors[i][2]*properties->shape_eigen_vectors[i][2]);
        for(j=0;j<3;j++)
          properties->shape_eigen_vectors[i][j]/=norm;
      }
    }
    // ... else apply defaults to faulty cases.
    else{ 
      if(flag_interpolated==FLAG_INTERP_MIN_BIN)      i_profile=0;
      else if(flag_interpolated==FLAG_INTERP_MAX_BIN) i_profile=profile->n_bins-1;
      else SID_trap_error("Unrecognized value for flag_interpolated {%d}.",flag_interpolated);

      //  ... COM positions ...
      properties->position_COM[0]=(double)(profile->bins[i_profile].position_COM[0])/expansion_factor;
      properties->position_COM[1]=(double)(profile->bins[i_profile].position_COM[1])/expansion_factor;
      properties->position_COM[2]=(double)(profile->bins[i_profile].position_COM[2])/expansion_factor;
      //  ... M_vir ...
      properties->M_vir=(double)profile->bins[i_profile].M_r;
      //  ... sigma_v ...
      properties->sigma_v=(float)profile->bins[i_profile].sigma_tot;
      //   ... spin ...
      properties->spin[0]=(float)profile->bins[i_profile].spin[0];
      properties->spin[1]=(float)profile->bins[i_profile].spin[1];
      properties->spin[2]=(float)profile->bins[i_profile].spin[2];
      //  ... triaxial axes ratios ...
      properties->q_triaxial=(float)profile->bins[i_profile].q_triaxial;
      properties->s_triaxial=(float)profile->bins[i_profile].s_triaxial;
      // ... shape eigen vectors ...
      for(i=0;i<3;i++){
        for(j=0;j<3;j++)
          properties->shape_eigen_vectors[i][j]=(float)profile->bins[i_profile].shape_eigen_vectors[i][j];
        norm=sqrt(properties->shape_eigen_vectors[i][0]*properties->shape_eigen_vectors[i][0]+
                  properties->shape_eigen_vectors[i][1]*properties->shape_eigen_vectors[i][1]+
                  properties->shape_eigen_vectors[i][2]*properties->shape_eigen_vectors[i][2]);
        for(j=0;j<3;j++)
          properties->shape_eigen_vectors[i][j]/=norm;
      }
    }

    // Enforce periodic box on COM position
    properties->position_COM[0]+=x_cen;
    properties->position_COM[1]+=y_cen;
    properties->position_COM[2]+=z_cen;
    if(properties->position_COM[0]< box_size) properties->position_COM[0]+=box_size;
    if(properties->position_COM[1]< box_size) properties->position_COM[1]+=box_size;
    if(properties->position_COM[2]< box_size) properties->position_COM[2]+=box_size;
    if(properties->position_COM[0]>=box_size) properties->position_COM[0]-=box_size;
    if(properties->position_COM[1]>=box_size) properties->position_COM[1]-=box_size;
    if(properties->position_COM[2]>=box_size) properties->position_COM[2]-=box_size;

    // Perform unit conversions
    //   ... properties first ...
    properties->position_COM[0]*=h_Hubble/M_PER_MPC;
    properties->position_COM[1]*=h_Hubble/M_PER_MPC;
    properties->position_COM[2]*=h_Hubble/M_PER_MPC;
    properties->position_MBP[0]*=h_Hubble/M_PER_MPC;
    properties->position_MBP[1]*=h_Hubble/M_PER_MPC;
    properties->position_MBP[2]*=h_Hubble/M_PER_MPC;
    properties->velocity_COM[0]*=1e-3;
    properties->velocity_COM[1]*=1e-3;
    properties->velocity_COM[2]*=1e-3;
    properties->velocity_MBP[0]*=1e-3;
    properties->velocity_MBP[1]*=1e-3;
    properties->velocity_MBP[2]*=1e-3;
    properties->M_vir          *=h_Hubble/M_SOL;
    properties->R_vir          *=h_Hubble/M_PER_MPC;
    properties->R_halo         *=h_Hubble/M_PER_MPC;
    properties->R_max          *=h_Hubble/M_PER_MPC;
    properties->V_max          *=1e-3;
    properties->sigma_v        *=1e-3;
    properties->spin[0]        *=1e-3*h_Hubble/M_PER_MPC;
    properties->spin[1]        *=1e-3*h_Hubble/M_PER_MPC;
    properties->spin[2]        *=1e-3*h_Hubble/M_PER_MPC;

    //   ... then profiles ...
    for(i_bin=0;i_bin<profile->n_bins;i_bin++){
      profile->bins[i_bin].r_med          *=h_Hubble/M_PER_MPC;
      profile->bins[i_bin].r_max          *=h_Hubble/M_PER_MPC;
      profile->bins[i_bin].M_r            *=h_Hubble/M_SOL;
      profile->bins[i_bin].rho            *=M_PER_MPC*M_PER_MPC*M_PER_MPC/(h_Hubble*h_Hubble*M_SOL);
      profile->bins[i_bin].position_COM[0]*=h_Hubble/M_PER_MPC;
      profile->bins[i_bin].position_COM[1]*=h_Hubble/M_PER_MPC;
      profile->bins[i_bin].position_COM[2]*=h_Hubble/M_PER_MPC;
      profile->bins[i_bin].velocity_COM[0]*=1e-3;
      profile->bins[i_bin].velocity_COM[1]*=1e-3;
      profile->bins[i_bin].velocity_COM[2]*=1e-3;
      profile->bins[i_bin].sigma_rad      *=1e-3;
      profile->bins[i_bin].sigma_tan      *=1e-3;
      profile->bins[i_bin].sigma_tot      *=1e-3;
      profile->bins[i_bin].spin[0]        *=1e-3*h_Hubble/M_PER_MPC;
      profile->bins[i_bin].spin[1]        *=1e-3*h_Hubble/M_PER_MPC;
      profile->bins[i_bin].spin[2]        *=1e-3*h_Hubble/M_PER_MPC;
    }
  
  }

  return(flag_interpolated);
}

