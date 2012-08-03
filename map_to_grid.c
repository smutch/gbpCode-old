#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpClustering.h>

void map_to_grid(size_t      n_particles_local, 
                 GBPREAL    *x_particles_local,
                 GBPREAL    *y_particles_local,
                 GBPREAL    *z_particles_local,
                 GBPREAL    *m_particles_local,
                 cosmo_info *cosmo,
                 double      redshift,
                 int         distribution_scheme,
                 double      normalization_target,
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
  int         flag_viable;
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
  int         W_search_lo;
  int         W_search_hi;
  size_t      receive_left_size=0;
  size_t      receive_right_size=0;
  size_t      index_best;
  field_info  buffer_left;
  field_info  buffer_right;
  int         n_buffer[3];
  size_t      n_send_left;
  size_t      n_send_right;
  size_t      send_size_left;
  size_t      send_size_right;
  GBPREAL    *send_left;
  GBPREAL    *send_right;
  GBPREAL    *receive_left=NULL;
  GBPREAL    *receive_right=NULL;
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

  // Compute the total poulation size and print a status message
  calc_sum_global(&n_particles_local,&n_particles,1,SID_SIZE_T,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  SID_log("Distributing %zu items onto a %dx%dx%d grid...",
          SID_LOG_OPEN,n_particles,field->n[0],field->n[1],field->n[2]);

  // Set some variables
  if(m_particles_local!=NULL)
    flag_multimass=TRUE;
  else{
    flag_multimass=FALSE;
    m_p=1.;
  }
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];

  // Initializing the mass assignment scheme
  switch(distribution_scheme){
  case MAP2GRID_DIST_DWT20:
    W_search_lo=2;
    W_search_hi=7;
    kernal_offset=2.5;
    compute_Daubechies_scaling_fctns(20,5,&r_Daub,&W_Daub,&n_Daub);
    init_interpolate(r_Daub,W_Daub,n_Daub,gsl_interp_cspline,&W_r_Daub_interp);
    SID_free(SID_FARG r_Daub);
    SID_free(SID_FARG W_Daub);
    SID_log("(using D20 scale function kernal)...",SID_LOG_CONTINUE);
    break;
  case MAP2GRID_DIST_DWT12:
    W_search_lo=1;
    W_search_hi=6;
    kernal_offset=1.75;
    compute_Daubechies_scaling_fctns(12,5,&r_Daub,&W_Daub,&n_Daub);
    init_interpolate(r_Daub,W_Daub,(size_t)n_Daub,gsl_interp_cspline,&W_r_Daub_interp);
    SID_free(SID_FARG r_Daub);
    SID_free(SID_FARG W_Daub);
    SID_log("(using D12 scale function kernal)...",SID_LOG_CONTINUE);
    break;
  case MAP2GRID_DIST_TSC:
    W_search_lo=2;
    W_search_hi=2;
    SID_log("(using triangular shaped function kernal)...",SID_LOG_CONTINUE);
    break;
  case MAP2GRID_DIST_CIC:
    SID_log("(using cloud-in-cell kernal)...",SID_LOG_CONTINUE);
  case MAP2GRID_DIST_NGP:
  default:
    W_search_lo=1;
    W_search_hi=1;
    SID_log("(using nearest grid point kernal)...",SID_LOG_CONTINUE);
    break;
  }

  // Initializing slab buffers
  n_send_left    =(size_t)(field->n[0]*field->n[1]*W_search_lo);
  n_send_right   =(size_t)(field->n[0]*field->n[1]*W_search_hi);
  send_size_left =n_send_left *sizeof(GBPREAL);
  send_size_right=n_send_right*sizeof(GBPREAL);
  send_left      =(GBPREAL *)SID_calloc(send_size_left);
  send_right     =(GBPREAL *)SID_calloc(send_size_right);
  receive_left   =(GBPREAL *)SID_calloc(send_size_right);
  receive_right  =(GBPREAL *)SID_calloc(send_size_left);

  // Create the mass distribution
  SID_log("Performing grid assignment...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Clear the field
  clear_field(field);

  // It is essential that we not pad the field for the simple way that we add-in the boundary buffers below
  set_FFT_padding_state(field,FALSE);

  // Loop over all the objects
  for(i_p=0,norm_local=0.,i_report=0,i_p_next_report=n_particles_local/10;i_p<n_particles_local;i_p++){     
    if(flag_multimass)
      m_p=(double)(m_particles_local[i_p]);

    // Particle's position
    x_particle_i=(GBPREAL)x_particles_local[i_p];
    y_particle_i=(GBPREAL)y_particles_local[i_p];
    z_particle_i=(GBPREAL)z_particles_local[i_p];

    // Quantize it onto the grid
    x_particle_i/=(GBPREAL)field->dR[0];
    y_particle_i/=(GBPREAL)field->dR[1];
    z_particle_i/=(GBPREAL)field->dR[2];
    i_i[0]=(int)x_particle_i; // position in grid-coordinates
    i_i[1]=(int)y_particle_i; // position in grid-coordinates
    i_i[2]=(int)z_particle_i; // position in grid-coordinates

    // Apply the kernel
    flag_viable=TRUE;
    double x_i_effective;
    for(j_i[0]=-W_search_lo;j_i[0]<=W_search_hi;j_i[0]++){
      for(j_i[1]=-W_search_lo;j_i[1]<=W_search_hi;j_i[1]++){
        for(j_i[2]=-W_search_lo;j_i[2]<=W_search_hi;j_i[2]++){
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
              x_i_effective=x_i+kernal_offset;
              if(x_i_effective>0.)
                W_i*=interpolate(W_r_Daub_interp,x_i_effective);
              else
                flag_active=FALSE;
              break;
              // Distribute using the triangular shaped cloud (TSC) method
            case MAP2GRID_DIST_TSC:
              if(x_i<0.5)
                W_i*=(0.75-x_i*x_i);
              else if(x_i<1.5)
                W_i*=0.5*(1.5-fabs(x_i))*(1.5-fabs(x_i));
              else
                flag_active=FALSE;
              break;
              // Distribute using the cloud-in-cell (CIC) method
            case MAP2GRID_DIST_CIC:
              if(fabs(x_i)<1.)
                W_i*=(1.-fabs(x_i));
              else
                flag_active=FALSE;
              break;
              // Distribute using "nearest grid point" (NGP; ie. the simplest and default) method
            case MAP2GRID_DIST_NGP:
            default:
              if(fabs(x_i)<=0.5 && flag_viable)
                W_i*=1.;
              else
                flag_active=FALSE;
              break;
            }
          }
          if(flag_active){ // This flags-out regions of the kernal with no support to save some time
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
            //     local array or to the slab buffers.
            k_i[0]=(i_i[0]+j_i[0]);
            if(k_i[0]<field->i_R_start_local[0]){
              k_i[0]-=(field->i_R_start_local[0]-W_search_lo);
              if(k_i[0]<0)
                 SID_trap_error("Left slab buffer limit exceeded by %d element(s).",ERROR_LOGIC,-k_i[0]);
              send_left[index_FFT_R(field,k_i)]+=W_i;
            }
            else if(k_i[0]>field->i_R_stop_local[0]){
              k_i[0]-=(field->i_R_stop_local[0]+1);
              if(k_i[0]>=W_search_hi)
                 SID_trap_error("Right slab buffer limit exceeded by %d element(s).",ERROR_LOGIC,k_i[0]-W_search_hi+1);
              send_right[index_FFT_R(field,k_i)]+=W_i;
            }
            else
              field->field_local[index_local_FFT_R(field,k_i)]+=W_i;
            norm_local+=W_i;
            flag_viable=FALSE;
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

  // Perform exchange of slab buffers and add them to the local mass distribution.
  //    Note: it's important that the FFT field not be padded (see above, where
  //          this is set) for this to work the way it's done.
  SID_log("Adding-in the slab buffers...",SID_LOG_OPEN|SID_LOG_TIMER);
  exchange_slab_buffer_left(send_left,
                            send_size_left,
                            receive_right,
                            &receive_right_size,
                            &(field->slab));
  exchange_slab_buffer_right(send_right,
                             send_size_right,
                             receive_left,
                             &receive_left_size,
                             &(field->slab));
  for(i_b=0;i_b<n_send_right;i_b++)
    field->field_local[i_b]+=receive_left[i_b];
  for(i_b=0;i_b<n_send_left;i_b++)
    field->field_local[field->n_field_R_local-n_send_left+i_b]+=receive_right[i_b];
  SID_free(SID_FARG send_left);
  SID_free(SID_FARG send_right);
  SID_free(SID_FARG receive_left);
  SID_free(SID_FARG receive_right);
  SID_log("Done.",SID_LOG_CLOSE);
  
  // Recompute local normalization (more accurate for large sample sizes)
  SID_log("Applying normalization...",SID_LOG_OPEN);
  norm_local=0;
  for(i_grid=0;i_grid<field->total_local_size;i_grid++)
    norm_local+=(double)field->field_local[i_grid];  
  calc_sum_global(&norm_local,&normalization,1,SID_DOUBLE,CALC_MODE_DEFAULT,SID.COMM_WORLD);
  double normalization_factor;
  normalization_factor=normalization/normalization_target;
  for(i_grid=0;i_grid<field->total_local_size;i_grid++)
    field->field_local[i_grid]*=normalization_factor;
  SID_log("Done.",SID_LOG_CLOSE,normalization);

  if(W_r_Daub_interp!=NULL)
    free_interpolate(SID_FARG W_r_Daub_interp);

  SID_log("Done.",SID_LOG_CLOSE);
  
}

