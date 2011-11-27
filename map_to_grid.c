#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>

double map_to_grid(size_t      n_particles_local, 
                   GBPREAL    *x_particles_local,
                   GBPREAL    *y_particles_local,
                   GBPREAL    *z_particles_local,
                   GBPREAL    *vx_particles_local,
                   GBPREAL    *vy_particles_local,
                   GBPREAL    *vz_particles_local,
                   GBPREAL    *m_particles_local,
                   cosmo_info *cosmo,
                   double      redshift,
                   int         distribution_scheme,
                   field_info *field){
  size_t      i_p;
  int         i_k;
  size_t      i_b;
  size_t      i_grid;
  int         i_coord;
  int         i_i[3];
  int         j_i[3];
  int         k_i[3];
  size_t      n_particles;
  double      m_p;
  int         flag_multimass;
  int         flag_active;
  int         flag_unused;
  double      k_mag;
  double      dk;
  int         n_powspec;
  int         mode_powspec;
  size_t     *n_mode_powspec;
  double     *k_powspec;
  double     *kmin_powspec;
  double     *kmax_powspec;
  double     *k_powspec_bin;
  double     *P_powspec;
  double     *dP_powspec;
  double      k_min;
  double      k_max;
  double      norm_local;
  double      normalization;
  double      x_i;
  GBPREAL     x_particle_i;
  GBPREAL     y_particle_i;
  GBPREAL     z_particle_i;
  double      kernal_offset;
  int         W_search;
  size_t      n_send;
  size_t      send_size;
  size_t      receive_left_size=0;
  size_t      receive_right_size=0;
  size_t      index_best;
  field_info  buffer_left;
  field_info  buffer_right;
  int         n_buffer[3];
  fftw_real  *send_left;
  fftw_real  *send_right;
  fftw_real  *receive_left=NULL;
  fftw_real  *receive_right=NULL;
  double       r_i,r_min,r_i_max=0;
  double       W_i;
  int          index_i;
  interp_info *P_k_interp;
  double      *r_Daub;
  double      *W_Daub;
  double       h_Hubble;
  int          n_Daub;
  size_t       i_p_next_report;
  int          i_report;
  interp_info *W_r_Daub_interp=NULL;
  int          i_rank;
  size_t       buffer_index;
  size_t      *field_sum_index;
  int          i_test;
  double       accumulator;

  calc_sum_global(&n_particles_local,&n_particles,1,SID_SIZE_T,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  SID_log("Distributing %zu items onto a %dx%dx%d grid...",
          SID_LOG_OPEN,n_particles,field->n[0],field->n[1],field->n[2]);

  if(m_particles_local!=NULL)
    flag_multimass=TRUE;
  else{
    flag_multimass=FALSE;
    m_p           =1.;
  }

  // Initializing the mass assignment scheme
  switch(distribution_scheme){
  case MAP2GRID_DIST_DWT20:
    W_search=4;
    kernal_offset=2.5;
    compute_Daubechies_scaling_fctns(20,5,&r_Daub,&W_Daub,&n_Daub);
    init_interpolate(r_Daub,W_Daub,n_Daub,gsl_interp_cspline,&W_r_Daub_interp);
    accumulator=0.;
    for(j_i[0]=-W_search+1;j_i[0]<=W_search;j_i[0]++){
      for(j_i[1]=-W_search+1;j_i[1]<=W_search;j_i[1]++){
        for(j_i[2]=-W_search+1;j_i[2]<=W_search;j_i[2]++){
          for(i_coord=0,W_i=1.;i_coord<3;i_coord++){
            switch(i_coord){
            case 0:
              x_i=(double)(j_i[0]);
              break;
            case 1:
              x_i=(double)(j_i[1]);
              break;
            case 2:
              x_i=(double)(j_i[2]);
              break;
            }
            if(fabs(x_i)<=(double)W_search)
              W_i*=interpolate(W_r_Daub_interp,x_i-4.);
          }
/*
          r_i=sqrt((double)(j_i[0]*j_i[0])+(double)(j_i[1]*j_i[1])+(double)(j_i[2]*j_i[2]));
          if(r_i<(double)W_search)
            W_i=interpolate(W_r_Daub_interp,r_i);
          else
            W_i=0.;
*/
          accumulator+=W_i;
        }
      }
    }
    fprintf(stderr,"test %le\n",accumulator);
    SID_free((void **)&r_Daub);
    SID_free((void **)&W_Daub);
    SID_log("(using D20 scale function kernal)...",SID_LOG_CONTINUE);
    break;
  case MAP2GRID_DIST_DWT12:
    W_search=3;
    kernal_offset=1.75;
    compute_Daubechies_scaling_fctns(12,5,&r_Daub,&W_Daub,&n_Daub);
    init_interpolate(r_Daub,W_Daub,(size_t)n_Daub,gsl_interp_cspline,&W_r_Daub_interp);
    SID_free((void **)&r_Daub);
    SID_free((void **)&W_Daub);
    SID_log("(using D12 scale function kernal)...",SID_LOG_CONTINUE);
    break;
  case MAP2GRID_DIST_TSC:
    W_search=2;
    SID_log("(using triangular shaped function kernal)...",SID_LOG_CONTINUE);
    break;
  case MAP2GRID_DIST_CIC:
    SID_log("(using cloud-in-cell kernal)...",SID_LOG_CONTINUE);
  case MAP2GRID_DIST_NGP:
  default:
    W_search=1;
    SID_log("(using nearest grid point kernal)...",SID_LOG_CONTINUE);
    break;
  }

  // Initializing slab buffers
  n_send       =(size_t)(field->n[0]*field->n[1]*W_search);
  send_size    =n_send*sizeof(fftw_real);
  send_left    =(fftw_real *)SID_malloc(send_size);
  send_right   =(fftw_real *)SID_malloc(send_size);
  receive_left =(fftw_real *)SID_malloc(send_size);
  receive_right=(fftw_real *)SID_malloc(send_size);
  for(i_b=0;i_b<n_send;i_b++){
    send_left[i_b] =0.;
    send_right[i_b]=0.;
  }

  // Create the mass distribution
  SID_log("Performing grid assignment...",SID_LOG_OPEN|SID_LOG_TIMER);
  clear_field(field);
  remove_buffer_FFT_R(field); // Essential for the simple way that we add-in the boundary buffers below
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  for(i_p=0,norm_local=0.,i_report=0,i_p_next_report=n_particles_local/10;i_p<n_particles_local;i_p++){     
    if(flag_multimass)
      m_p=(double)(m_particles_local[i_p]);
    x_particle_i=(GBPREAL)x_particles_local[i_p];
    y_particle_i=(GBPREAL)y_particles_local[i_p];
    z_particle_i=(GBPREAL)z_particles_local[i_p];
    if(vx_particles_local!=NULL)
      x_particle_i+=(GBPREAL)(1e3*h_Hubble*((double)vx_particles_local[i_p])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
    if(vy_particles_local!=NULL)
      y_particle_i+=(GBPREAL)(1e3*h_Hubble*((double)vy_particles_local[i_p])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
    if(vz_particles_local!=NULL)
      z_particle_i+=(GBPREAL)(1e3*h_Hubble*((double)vz_particles_local[i_p])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
    force_periodic(&x_particle_i,0.,field->L[0]);
    force_periodic(&y_particle_i,0.,field->L[1]);
    force_periodic(&z_particle_i,0.,field->L[2]);
    x_particle_i/=(GBPREAL)field->dR[0];
    y_particle_i/=(GBPREAL)field->dR[1];
    z_particle_i/=(GBPREAL)field->dR[2];
    i_i[0]=(int)x_particle_i; // position in grid-coordinates
    i_i[1]=(int)y_particle_i; // position in grid-coordinates
    i_i[2]=(int)z_particle_i; // position in grid-coordinates
    flag_unused=TRUE;
    for(j_i[0]=-W_search+1;j_i[0]<=W_search;j_i[0]++){
      for(j_i[1]=-W_search+1;j_i[1]<=W_search;j_i[1]++){
        for(j_i[2]=-W_search+1;j_i[2]<=W_search;j_i[2]++){
          // Compute distance to each grid point being searched against ...
          flag_active=TRUE;
          for(i_coord=0,W_i=1.;i_coord<3;i_coord++){
            switch(i_coord){
            case 0:
              x_i=(double)(i_i[0]+j_i[0])-(double)x_particle_i;
              break;
            case 1:
              x_i=(double)(i_i[1]+j_i[1])-(double)y_particle_i;
              break;
            case 2:
              x_i=(double)(i_i[2]+j_i[2])-(double)z_particle_i;
              break;
            }
            switch(distribution_scheme){
              // Distribute with a Daubechies wavelet transform of 12th or 20th order a la Cui et al '08
            case MAP2GRID_DIST_DWT12:
            case MAP2GRID_DIST_DWT20:
              if(fabs(x_i)<=(double)W_search)
                W_i*=interpolate(W_r_Daub_interp,x_i);
              else{
                W_i=0.;
                flag_active=FALSE;
              }
              break;
              // Distribute using the triangular shaped cloud (TSC) method
            case MAP2GRID_DIST_TSC:
              if(x_i<0.5)
                W_i*=(0.75-x_i*x_i);
              else if(x_i<1.5)
                W_i*=0.5*(1.5-fabs(x_i))*(1.5-fabs(x_i));
              else{
                W_i=0.;
                flag_active=FALSE;
              }
              break;
              // Distribute using the cloud-in-cell (CIC) method
            case MAP2GRID_DIST_CIC:
              if(fabs(x_i)<1.)
                W_i*=(1.-fabs(x_i));
              else{
                W_i=0.;
                flag_active=FALSE;
              }
              break;
              // Distribute using "nearest grid point" (NGP; ie. the simplest and default) method
            case MAP2GRID_DIST_NGP:
            default:
              if(fabs(x_i)<=0.5 && flag_unused)
                W_i*=1.;
              else{
                W_i=0.;
                flag_active=FALSE;
              }
              break;
            }
          }
          if(flag_active){ // W_i can be negative, so we use this flag to decide when we don't need to do this
            W_i*=m_p;
            // Set the grid indices (enforce periodic BCs; do x-coordinate last) ...
            //   ... y-coordinate ...
            k_i[1]=(i_i[1]+j_i[1]);
            if(k_i[1]<0)
              k_i[1]+=field->n[1];
            else
              k_i[1]=k_i[1]%field->n[1];
            //   ... z-coordinate ...
            k_i[2]=i_i[2]+j_i[2];
            if(k_i[2]<0)
              k_i[2]+=field->n[2];
            else 
              k_i[2]=k_i[2]%field->n[2];
            //   ... x-coordinate ... 
            //     Depending on x-index, add contribution to the
            //     local array or to the slab buffers
            k_i[0]=(i_i[0]+j_i[0]);
            k_i[0]-=field->i_R_start_local[0];
            if(k_i[0]<0){
              k_i[0]+=W_search;
              send_left[index_FFT_R(field,k_i)]+=W_i;
            }
            else if(k_i[0]>=field->n_R_local[0]){
              k_i[0]-=field->n_R_local[0];
              send_right[index_FFT_R(field,k_i)]+=W_i;
            }
            else
              field->field_local[index_FFT_R(field,k_i)]+=W_i;
            norm_local+=W_i;
            flag_unused=FALSE;
          }
        }
      }
    }
    if(SID.I_am_Master){
     if(i_p==i_p_next_report){
       i_report++;
       SID_log("%3d%% complete.",SID_LOG_COMMENT|SID_LOG_TIMER,10*i_report);
       i_p_next_report=MIN(n_particles_local,n_particles_local*(i_report+1)/10);
     } 
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Perform exchange of slab buffers and add them to the local mass distribution
  SID_log("Adding-in the slab buffers...",SID_LOG_OPEN|SID_LOG_TIMER);
  exchange_slab_buffer_left(send_left,
                            send_size,
                            receive_right,
                            &receive_right_size,
                            &(field->slab));
  exchange_slab_buffer_right(send_right,
                             send_size,
                             receive_left,
                             &receive_left_size,
                             &(field->slab));
  for(i_b=0;i_b<n_send;i_b++)
    field->field_local[i_b]+=receive_left[i_b];
  for(i_b=0;i_b<n_send;i_b++)
    field->field_local[field->n_field_R_local-n_send+i_b]+=receive_right[i_b];
  SID_free((void **)&send_left);
  SID_free((void **)&send_right);
  SID_free((void **)&receive_left);
  SID_free((void **)&receive_right);
  SID_log("Done.",SID_LOG_CLOSE);
  
  // Recompute local normalization (more accurate for large sample sizes)
  SID_log("Computing normalization...",SID_LOG_OPEN);
  norm_local=0;
  for(i_grid=0;i_grid<field->total_local_size;i_grid++)
    norm_local+=(double)field->field_local[i_grid];  
  calc_sum_global(&norm_local,&normalization,1,SID_DOUBLE,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  SID_log("Done. (normalization=%le)",SID_LOG_CLOSE,normalization);

  if(W_r_Daub_interp!=NULL)
    free_interpolate(SID_FARG W_r_Daub_interp);

  SID_log("Done.",SID_LOG_CLOSE);
  
  return(normalization);
}

