#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpClustering.h>

void compute_pspec(plist_info  *plist,
                   char        *species_name,
                   pspec_info  *pspec,
                   int          i_run){
  int         i_k;
  int         i_coord;
  int         i_i[3];
  int         i_i_x[3];
  int         i_i_y[3];
  int         j_i[3];
  int         k_i[3];
  char        N_name[256];
  char        n_name[256];
  char        x_name[256];
  char        y_name[256];
  char        z_name[256];
  char        m_name[256];
  char        marray_name[256];
  char        i_name[256];
  char        k_name[256];
  char        P_name[256];
  char        store_name[256];
  char        coord_name[20];
  size_t      n_particles_local;
  size_t      n_particles;
  size_t     *cell_index_local;
  GBPREAL       *x_particles_local;
  GBPREAL       *y_particles_local;
  GBPREAL       *z_particles_local;
  GBPREAL       *vx_particles_local;
  GBPREAL       *vy_particles_local;
  GBPREAL       *vz_particles_local;
  GBPREAL       *m_particles_local;
  double      m_p;
  double     *k_2D;
  double     *k_2D_local;
  int         flag_multimass;
  int         flag_active;
  double      k_mag;
  double      dk;
  double      k_mag_x,k_mag_y;
  int         n_k;
  int         mode_powspec;
  int         mode_x,mode_y;
  int        *n_modes_1D_local;
  int        *n_modes_2D_local;
  double     *k_1D_local;
  double     *P_k_1D_local;
  double     *P_k_2D_local;
  double     *sigma_P_powspec;
  double     *sigma_P_powspec_local;
  double     *dP_powspec;
  double      shot_noise;
  double      k_min;
  double      k_max;
  double      norm_local;
  double      normalization;
  double      x_i;
  int         W_search;
  size_t      n_send;
  size_t      send_size;
  size_t      receive_left_size=0;
  size_t      receive_right_size=0;
  size_t      index_best;
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
  double       shot_noise_arg;
  int          n_Daub;
  interp_info *W_r_Daub_interp;
  FILE        *fp_test;
  double       dk_1D;
  double       dk_2D;

  SID_log("Computing power spectrum...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Parse some stuff from the pspec structure
  field_info  *FFT;
  int          distribution_scheme;
  cosmo_info  *cosmo;
  double       redshift;
  int          n_k_1D;
  double       k_min_1D;
  double       k_max_1D;
  int          n_k_2D;
  double       k_min_2D;
  double       k_max_2D;
  double      *k_1D;
  double      *P_k_1D;
  double      *dP_k_1D;
  int         *n_modes_1D;
  double      *P_k_2D;
  double      *dP_k_2D;
  int         *n_modes_2D;
  FFT                =&(pspec->FFT);
  distribution_scheme=pspec->mass_assignment_scheme;
  cosmo              =pspec->cosmo;
  redshift           =pspec->redshift;
  n_k_1D             =pspec->n_k_1D;
  k_min_1D           =pspec->k_min_1D;
  k_max_1D           =pspec->k_max_1D;
  n_k_2D             =pspec->n_k_2D;
  k_min_2D           =pspec->k_min_2D;
  k_max_2D           =pspec->k_max_2D;
  k_1D               =pspec->k_1D;
  P_k_1D             =pspec->P_k_1D[i_run];
  dP_k_1D            =pspec->dP_k_1D[i_run];
  n_modes_1D         =pspec->n_modes_1D;
  P_k_2D             =pspec->P_k_2D[i_run];
  dP_k_2D            =pspec->dP_k_2D[i_run];
  n_modes_2D         =pspec->n_modes_2D;

  // Decide if we are computing real or redshift space results
  int mode;
  switch(i_run){
  case 0:
     mode=PSPEC_DEFAULT;
     break;
  case 1:
     mode=PSPEC_ADD_VX;
     break;
  case 2:
     mode=PSPEC_ADD_VY;
     break;
  case 3:
     mode=PSPEC_ADD_VZ;
     break;
  }

  dk_1D=(k_max_1D-k_min_1D)/(double)(n_k_1D);
  dk_2D=(k_max_2D-k_min_2D)/(double)(n_k_2D);

  // Fetch the needed information
  n_particles        =((size_t    *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0]; 
  n_particles_local  =((size_t    *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0]; 
  x_particles_local  = (GBPREAL   *)ADaPS_fetch(plist->data,"x_%s",    species_name);
  y_particles_local  = (GBPREAL   *)ADaPS_fetch(plist->data,"y_%s",    species_name);
  z_particles_local  = (GBPREAL   *)ADaPS_fetch(plist->data,"z_%s",    species_name);
  if(ADaPS_exist(plist->data,"M_%s",species_name))
    m_particles_local=(GBPREAL *)ADaPS_fetch(plist->data,"M_%s",species_name);
  else
    m_particles_local=NULL;

  // Generate mass-field
  map_to_grid(n_particles_local, 
              x_particles_local,
              y_particles_local,
              z_particles_local,
              m_particles_local,
              cosmo,
              redshift,
              distribution_scheme,
              (double)n_particles,
              FFT);

  // Compute the FFT of the mass-field
  SID_log("Computing FFT...",SID_LOG_OPEN|SID_LOG_TIMER);
  compute_FFT(FFT);
  SID_log("Done.",SID_LOG_CLOSE);

  // Allocate local arrays; Initialize them and global arrays where results are stores
  k_1D_local      =(double *)SID_calloc(sizeof(double)*(n_k_1D)); 
  P_k_1D_local    =(double *)SID_calloc(sizeof(double)*(n_k_1D)); 
  n_modes_1D_local=(int    *)SID_calloc(sizeof(int)*(n_k_1D)); 
  for(i_k=0;i_k<n_k_1D;i_k++){
    P_k_1D[i_k]    =0.;
    dP_k_1D[i_k]   =0.;
    n_modes_1D[i_k]=0;
  }

  P_k_2D_local    =(double *)SID_calloc(sizeof(double)*(n_k_2D)*(n_k_2D));
  n_modes_2D_local=(int    *)SID_calloc(sizeof(int)*(n_k_2D)*(n_k_2D));
  for(i_k=0;i_k<n_k_2D*n_k_2D;i_k++){
    P_k_2D[i_k]    =0.;
    dP_k_2D[i_k]   =0.;
    n_modes_2D[i_k]=0;
  }

  // Perform binning for 1D power spectrum
  SID_log("Performing 1-d bining of power spectrum...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_i[0]=FFT->i_k_start_local[0],j_i[0]=0;i_i[0]<=FFT->i_k_stop_local[0];i_i[0]++,j_i[0]++){
    for(i_i[1]=FFT->i_k_start_local[1],j_i[1]=0;i_i[1]<=FFT->i_k_stop_local[1];i_i[1]++,j_i[1]++){
      for(i_i[2]=FFT->i_k_start_local[2],j_i[2]=0;i_i[2]<=FFT->i_k_stop_local[2];i_i[2]++,j_i[2]++){
        k_mag       =k_mag_field_FFT(FFT,i_i);
        mode_powspec=(int)((k_mag-k_min_1D)/dk_1D);
        if(mode_powspec<(n_k_1D) && mode_powspec>=0){
          k_1D_local[mode_powspec]  +=k_mag;
          P_k_1D_local[mode_powspec]+=(pow((double)FFT->cfield_local[index_FFT_k(FFT,j_i)].re,2.)+
                                       pow((double)FFT->cfield_local[index_FFT_k(FFT,j_i)].im,2.));
          n_modes_1D_local[mode_powspec]++;
        }
      }
    }
  }

  // Normalize the quantities we are averaging
  for(i_k=0;i_k<n_k_1D;i_k++){
    calc_sum_global(&(n_modes_1D_local[i_k]),&(n_modes_1D[i_k]),1,SID_INT,   CALC_MODE_DEFAULT,SID.COMM_WORLD);
    calc_sum_global(&(k_1D_local[i_k]),      &(k_1D[i_k]),      1,SID_DOUBLE,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    calc_sum_global(&(P_k_1D_local[i_k]),    &(P_k_1D[i_k]),    1,SID_DOUBLE,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    k_1D[i_k]  /=(double)((n_modes_1D)[i_k]);
    P_k_1D[i_k]/=(double)((n_modes_1D)[i_k]);
  }

  // Deal with shot noise.  This differs, depending on the
  //   mass-assignment scheme we used (see Cui et al 2008).
  for(i_k=0;i_k<(n_k_1D);i_k++){
    P_k_1D[i_k] *=FFT->L[0]*FFT->L[1]*FFT->L[2]/pow((double)n_particles,2.);
    dP_k_1D[i_k] =(P_k_1D)[i_k]/sqrt(n_modes_1D[i_k]);
    switch(distribution_scheme){
    case MAP2GRID_DIST_CIC:
      shot_noise_arg=sin(M_PI*(k_1D)[i_k]/(2.*FFT->k_Nyquist[0]));
      shot_noise    =(1.-2.*pow(shot_noise_arg,2.)/3.)/(double)n_particles;
      break;
    case MAP2GRID_DIST_TSC:
      shot_noise_arg=sin(M_PI*(k_1D)[i_k]/(2.*FFT->k_Nyquist[0]));
      shot_noise    =(1.-pow(shot_noise_arg,2.)+2.*pow(shot_noise_arg,4.)/15.)/(double)n_particles;
      break;
    case MAP2GRID_DIST_DWT12:
    case MAP2GRID_DIST_DWT20:
    case MAP2GRID_DIST_NGP:
    default:
      shot_noise=1./(double)n_particles;
    }
    shot_noise  *=FFT->L[0]*FFT->L[1]*FFT->L[2];
    P_k_1D[i_k] -=shot_noise;
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Perform binning for 2D power spectrum
  SID_log("Performing 2-d bining of power spectrum...",SID_LOG_OPEN|SID_LOG_TIMER);
  k_2D_local=(double *)SID_malloc(sizeof(double)*((n_k_2D)*(n_k_2D)));
  k_2D      =(double *)SID_malloc(sizeof(double)*((n_k_2D)*(n_k_2D)));
  for(i_k=0;i_k<(n_k_2D)*(n_k_2D);i_k++){
    k_2D_local[i_k]=0.;
    k_2D[i_k]      =0.;
  }
  for(i_i[0]=FFT->i_k_start_local[0],j_i[0]=0;i_i[0]<=FFT->i_k_stop_local[0];i_i[0]++,j_i[0]++){
    for(i_i[1]=FFT->i_k_start_local[1],j_i[1]=0;i_i[1]<=FFT->i_k_stop_local[1];i_i[1]++,j_i[1]++){
      for(i_i[2]=FFT->i_k_start_local[2],j_i[2]=0;i_i[2]<=FFT->i_k_stop_local[2];i_i[2]++,j_i[2]++){
        if(check_mode_for_flag(mode,PSPEC_ADD_VX)){
          k_mag_x=sqrt(FFT->k_field[1][i_i[1]]*FFT->k_field[1][i_i[1]]+FFT->k_field[2][i_i[2]]*FFT->k_field[2][i_i[2]]);
          k_mag_y=FFT->k_field[0][i_i[0]];
        }
        else if(check_mode_for_flag(mode,PSPEC_ADD_VY)){
          k_mag_x=sqrt(FFT->k_field[0][i_i[0]]*FFT->k_field[0][i_i[0]]+FFT->k_field[2][i_i[2]]*FFT->k_field[2][i_i[2]]);
          k_mag_y=FFT->k_field[1][i_i[1]];
        }
        else{
          k_mag_x=sqrt(FFT->k_field[0][i_i[0]]*FFT->k_field[0][i_i[0]]+FFT->k_field[1][i_i[1]]*FFT->k_field[1][i_i[1]]);
          k_mag_y=FFT->k_field[2][i_i[2]];
        }
        mode_x =(int)((k_mag_x-k_min_2D)/dk_2D);
        mode_y =(int)((k_mag_y-k_min_2D)/dk_2D);
        if(mode_x<(n_k_2D) && mode_x>=0 && mode_y<(n_k_2D) && mode_y>=0){
          mode_powspec=mode_y*(n_k_2D)+mode_x;
          k_2D_local[mode_powspec]  +=k_mag_field_FFT(FFT,i_i);
          P_k_2D_local[mode_powspec]+=(pow((double)FFT->cfield_local[index_FFT_k(FFT,j_i)].re,2.)+
                                       pow((double)FFT->cfield_local[index_FFT_k(FFT,j_i)].im,2.));
          n_modes_2D_local[mode_powspec]++;
        }
      }
    }
  }
  for(i_k=0;i_k<(n_k_2D)*(n_k_2D);i_k++){
    calc_sum_global(&(n_modes_2D_local[i_k]),&(n_modes_2D[i_k]),1,SID_INT,   CALC_MODE_DEFAULT,SID.COMM_WORLD);
    calc_sum_global(&(k_2D_local[i_k]),      &(k_2D[i_k]),      1,SID_DOUBLE,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    calc_sum_global(&(P_k_2D_local[i_k]),    &(P_k_2D[i_k]),    1,SID_DOUBLE,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    k_2D[i_k]  /=(double)((n_modes_2D)[i_k]);
    P_k_2D[i_k]/=(double)((n_modes_2D)[i_k]);
  }
  for(i_k=0;i_k<(n_k_2D)*(n_k_2D);i_k++){
    P_k_2D[i_k] *=FFT->L[0]*FFT->L[1]*FFT->L[2]/pow(n_particles,2.);
    dP_k_2D[i_k] =(P_k_2D)[i_k]/sqrt(n_modes_2D[i_k]);
    switch(distribution_scheme){
    case MAP2GRID_DIST_CIC:
      shot_noise_arg=sin(M_PI*k_2D[i_k]/(2.*FFT->k_Nyquist[0]));
      shot_noise    =(1.-2.*pow(shot_noise_arg,2.)/3.)/(double)n_particles;
      break;
    case MAP2GRID_DIST_TSC:
      shot_noise_arg=sin(M_PI*k_2D[i_k]/(2.*FFT->k_Nyquist[0]));
      shot_noise    =(1.-pow(shot_noise_arg,2.)+2.*pow(shot_noise_arg,4.)/15.)/(double)n_particles;
      break;
    case MAP2GRID_DIST_DWT12:
    case MAP2GRID_DIST_DWT20:
    case MAP2GRID_DIST_NGP:
    default:
      shot_noise=1./(double)n_particles;
    }
    shot_noise *=FFT->L[0]*FFT->L[1]*FFT->L[2];
    P_k_2D[i_k] -=shot_noise;
  }
  SID_log("Done.",SID_LOG_CLOSE);
  
  SID_free(SID_FARG k_1D_local);
  SID_free(SID_FARG k_2D);
  SID_free(SID_FARG k_2D_local);
  SID_free(SID_FARG P_k_1D_local);
  SID_free(SID_FARG n_modes_1D_local);
  SID_free(SID_FARG P_k_2D_local);
  SID_free(SID_FARG n_modes_2D_local);

  // Tell the datastructure that the calculation is done
  pspec->flag_processed[i_run]=TRUE;

  SID_log("Done.",SID_LOG_CLOSE);

}

