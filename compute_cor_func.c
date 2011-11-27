#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpClustering.h>

double calc_CFUNC_local(double DD,double DR,double RR);
double calc_CFUNC_local(double DD,double DR,double RR){
  if(RR>0)
    return((DD-2*DR+RR)/RR);
  else
    return(0.);
}

void compute_cor_func(plist_info  *plist,
                      int          mode,
                      cosmo_info  *cosmo,
                      char        *species_name,
                      char        *random_name,
                      double       redshift,
                      double       box_size,
                      double       r_min_l1D,
                      double       r_max_1D,
                      double       r_min_2D,
                      double       r_max_2D,
                      int          n_1D,
                      int          n_2D,
                      int          n_jack,
                      double      *CFUNC_l1D,
                      double      *dCFUNC_l1D,
                      double      *COVMTX_l1D,
                      double      *CFUNC_1D,
                      double      *dCFUNC_1D,
                      double      *COVMTX_1D,
                      double      *CFUNC_2D,
                      double      *dCFUNC_2D,
                      double      *COVMTX_2D,
                      int         *flag_compute_RR,
                      long long  **DD_l1D,
                      long long  **DR_l1D,
                      long long  **RR_l1D,
                      long long  **DD_1D,
                      long long  **DR_1D,
                      long long  **RR_1D,
                      long long  **DD_2D,
                      long long  **DR_2D,
                      long long  **RR_2D){
  int         i_bin,j_bin;
  int         i_jack;
  int         i_x_jack,i_y_jack,i_z_jack;
  int         n_jack_total;
  int         n_2D_total;
  int         i_rank;
  int         j_rank;
  size_t      i_data;
  size_t      j_data;
  size_t      n_data_local;
  size_t      n_data;
  size_t      i_random;
  size_t      j_random;
  size_t      n_random_local;
  size_t      n_random;
  size_t      n_data_allocate;
  size_t      n_random_allocate;
  size_t      n_data_rank;
  size_t      n_temp;
  GBPREAL       *x_data_rank;
  GBPREAL       *y_data_rank;
  GBPREAL       *z_data_rank;
  size_t      n_random_rank;
  GBPREAL       *x_random_rank;
  GBPREAL       *y_random_rank;
  GBPREAL       *z_random_rank;
  size_t     *x_data_rank_index;
  size_t     *x_random_rank_index;
  GBPREAL       *x_random_local_swap;
  GBPREAL       *y_random_local_swap;
  GBPREAL       *z_random_local_swap;
  GBPREAL       *x_data_local_swap;
  GBPREAL       *y_data_local_swap;
  GBPREAL       *z_data_local_swap;
  GBPREAL       *x_data_local_temp;
  GBPREAL       *y_data_local_temp;
  GBPREAL       *z_data_local_temp;
  GBPREAL       *x_data_local;
  GBPREAL       *y_data_local;
  GBPREAL       *z_data_local;
  GBPREAL       *vx_data_local;
  GBPREAL       *vy_data_local;
  GBPREAL       *vz_data_local;
  GBPREAL       *x_random_local;
  GBPREAL       *y_random_local;
  GBPREAL       *z_random_local;
  GBPREAL        x_i,y_i,z_i;
  GBPREAL        x_j,y_j,z_j;
  int         flag_alloc_x;
  int         flag_alloc_y;
  int         flag_alloc_z;
  int         flag_alloc_x_swap;
  int         flag_alloc_y_swap;
  int         flag_alloc_z_swap;
  double      mean,std_dev;
  double      dx,dy,dz;
  double      sep_1D;
  double      sep_2D_x;
  double      sep_2D_y;
  int         bin_1D;
  int         bin_l1D;
  int         bin_2D_x;
  int         bin_2D_y;
  double      d_jack;
  int        *data_zone_local;
  int        *random_zone_local;
  int        *data_zone_rank;
  int        *random_zone_rank;
  size_t     *x_data_local_index;
  size_t     *x_random_local_index;
  size_t      j_data_lo_1;
  size_t      j_data_hi_1;
  size_t      j_data_lo_2;
  size_t      j_data_hi_2;
  size_t      j_random_lo_1;
  size_t      j_random_hi_1;
  size_t      j_random_lo_2;
  size_t      j_random_hi_2;
  int         bin_2D;
  double      r_max;
  double      r_min_1D;
  double      r_max_l1D;
  double      dr_l1D;
  double      dr_1D;
  double      dr_2D;
  GBPREAL        x_temp;
  GBPREAL        y_temp;
  GBPREAL        z_temp;
  long long  *temp_array;
  double      bar_i;
  double      bar_j;
  double      bar_ij;
  double      h_Hubble;

  SID_trap_error("This code is broken.  Pairs are being counted multiple times.  There is a problem with the last JK region.",ERROR_LOGIC);

  SID_log("Computing correlation function (%d 1D bins and %d 2D bins)...",SID_LOG_OPEN|SID_LOG_TIMER,n_1D,n_2D);
  r_max    =MAX(r_max_1D,r_max_2D);
  r_min_1D =0.;
  r_min_l1D=take_log10(r_min_l1D);
  r_max_l1D=take_log10(r_max_1D);
  dr_l1D   =(r_max_l1D-r_min_l1D)/(double)n_1D;
  dr_1D    =(r_max_1D)/(double)n_1D;
  dr_2D    =(r_max_2D-r_min_2D)/(double)n_2D;
  n_jack_total=n_jack*n_jack*n_jack;
  n_2D_total  =n_2D*n_2D;

  // Fetch the needed information
  SID_log("Initializing...",SID_LOG_OPEN);
  n_data           =((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0]; 
  n_data_local     =((size_t *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0]; 
  x_data_local_temp= (GBPREAL   *)ADaPS_fetch(plist->data,"x_%s",    species_name);
  y_data_local_temp= (GBPREAL   *)ADaPS_fetch(plist->data,"y_%s",    species_name);
  z_data_local_temp= (GBPREAL   *)ADaPS_fetch(plist->data,"z_%s",    species_name);
  flag_alloc_x     =FALSE;
  flag_alloc_y     =FALSE;
  flag_alloc_z     =FALSE;
  h_Hubble=((double *)ADaPS_fetch(cosmo,"h_Hubble"))[0];
  if(check_mode_for_flag(mode,CFUNC_ADD_VX)){
    vx_data_local=(GBPREAL *)ADaPS_fetch(plist->data,"vx_%s",species_name);
    x_data_local =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_data_local);
    for(i_data=0;i_data<n_data_local;i_data++){
      x_temp=x_data_local_temp[i_data]+(GBPREAL)(1e3*h_Hubble*((double)vx_data_local[i_data])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
      force_periodic(&x_temp,0.,(GBPREAL)box_size);
      x_data_local[i_data]=(GBPREAL)x_temp;
    }
    flag_alloc_x =TRUE;
  }
  else
    x_data_local=x_data_local_temp;
  if(check_mode_for_flag(mode,CFUNC_ADD_VY)){
    vy_data_local=(GBPREAL *)ADaPS_fetch(plist->data,"vy_%s",species_name);
    y_data_local =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_data_local);
    for(i_data=0;i_data<n_data_local;i_data++){
      y_temp=y_data_local_temp[i_data]+(GBPREAL)(1e3*h_Hubble*((double)vy_data_local[i_data])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
      force_periodic(&y_temp,0.,(GBPREAL)box_size);
      y_data_local[i_data]=(GBPREAL)y_temp;
    }
    flag_alloc_y =TRUE;
  }
  else
    y_data_local=y_data_local_temp;
  if(check_mode_for_flag(mode,CFUNC_ADD_VZ)){
    vz_data_local=(GBPREAL *)ADaPS_fetch(plist->data,"vz_%s",species_name);
    z_data_local =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_data_local);
    for(i_data=0;i_data<n_data_local;i_data++){
      z_temp=z_data_local_temp[i_data]+(GBPREAL)(1e3*h_Hubble*((double)vz_data_local[i_data])/(a_of_z(redshift)*M_PER_MPC*H_convert(H_z(redshift,cosmo))));
      force_periodic(&z_temp,0.,(GBPREAL)box_size);
      z_data_local[i_data]=(GBPREAL)z_temp;
    }
    flag_alloc_z =TRUE;
  }
  else
    z_data_local=z_data_local_temp;

  n_random      =((size_t *)ADaPS_fetch(plist->data,"n_all_%s",random_name))[0]; 
  n_random_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",    random_name))[0]; 
  x_random_local= (GBPREAL   *)ADaPS_fetch(plist->data,"x_%s",    random_name);
  y_random_local= (GBPREAL   *)ADaPS_fetch(plist->data,"y_%s",    random_name);
  z_random_local= (GBPREAL   *)ADaPS_fetch(plist->data,"z_%s",    random_name);

  // Make the z-coordinate the redshift-space coordinate
  flag_alloc_x_swap  =flag_alloc_x;
  flag_alloc_y_swap  =flag_alloc_y;
  flag_alloc_z_swap  =flag_alloc_z;
  x_data_local_swap  =x_data_local;
  y_data_local_swap  =y_data_local;
  z_data_local_swap  =z_data_local;
  x_random_local_swap=x_random_local;
  y_random_local_swap=y_random_local;
  z_random_local_swap=z_random_local;
  if(check_mode_for_flag(mode,CFUNC_ADD_VX)){
    flag_alloc_x  =flag_alloc_y_swap;
    flag_alloc_y  =flag_alloc_z_swap;
    flag_alloc_z  =flag_alloc_x_swap;
    x_data_local  =y_data_local_swap;
    y_data_local  =z_data_local_swap;
    z_data_local  =x_data_local_swap;
    x_random_local=y_random_local_swap;
    y_random_local=z_random_local_swap;
    z_random_local=x_random_local_swap;
  }
  else if(check_mode_for_flag(mode,CFUNC_ADD_VY)){
    flag_alloc_x  =flag_alloc_x_swap;
    flag_alloc_y  =flag_alloc_z_swap;
    flag_alloc_z  =flag_alloc_y_swap;
    x_data_local  =x_data_local_swap;
    y_data_local  =z_data_local_swap;
    z_data_local  =y_data_local_swap;
    x_random_local=x_random_local_swap;
    y_random_local=z_random_local_swap;
    z_random_local=y_random_local_swap;
  }
  else{
    flag_alloc_x  =flag_alloc_x_swap;
    flag_alloc_y  =flag_alloc_y_swap;
    flag_alloc_z  =flag_alloc_z_swap;
    x_data_local  =x_data_local_swap;
    y_data_local  =y_data_local_swap;
    z_data_local  =z_data_local_swap;
    x_random_local=x_random_local_swap;
    y_random_local=y_random_local_swap;
    z_random_local=z_random_local_swap;    
  }

  // Determine which jacknife zone each object belongs to
  d_jack=box_size/(double)n_jack;
  data_zone_local=(int *)SID_malloc(sizeof(int)*n_data_local);
  for(i_data=0;i_data<n_data_local;i_data++){
    i_x_jack=(int)(x_data_local[i_data]/d_jack);
    i_y_jack=(int)(y_data_local[i_data]/d_jack);
    i_z_jack=(int)(z_data_local[i_data]/d_jack);
    if(i_x_jack<0)       i_x_jack=0;
    if(i_x_jack>=n_jack) i_x_jack=n_jack-1;
    if(i_y_jack<0)       i_y_jack=0;
    if(i_y_jack>=n_jack) i_y_jack=n_jack-1;
    if(i_z_jack<0)       i_z_jack=0;
    if(i_z_jack>=n_jack) i_z_jack=n_jack-1;
    data_zone_local[i_data]=(i_x_jack*n_jack*n_jack)+i_y_jack*n_jack+i_z_jack;
  }
  random_zone_local=(int *)SID_malloc(sizeof(int)*n_random_local);
  for(i_random=0;i_random<n_random_local;i_random++){
    i_x_jack=(int)(x_random_local[i_random]/d_jack);
    i_y_jack=(int)(y_random_local[i_random]/d_jack);
    i_z_jack=(int)(z_random_local[i_random]/d_jack);
    if(i_x_jack<0)       i_x_jack=0;
    if(i_x_jack>=n_jack) i_x_jack=n_jack-1;
    if(i_y_jack<0)       i_y_jack=0;
    if(i_y_jack>=n_jack) i_y_jack=n_jack-1;
    if(i_z_jack<0)       i_z_jack=0;
    if(i_z_jack>=n_jack) i_z_jack=n_jack-1;
    random_zone_local[i_random]=(i_x_jack*n_jack*n_jack)+i_y_jack*n_jack+i_z_jack;
  }

  // Allocate and initialize temporary arrays for pair counts
  if((*flag_compute_RR)){
    for(i_jack=0;i_jack<=n_jack_total;i_jack++){
      for(i_bin=0;i_bin<n_1D;i_bin++){
        RR_1D[i_jack][i_bin] =0;
        RR_l1D[i_jack][i_bin]=0;
      }
      for(i_bin=0;i_bin<n_2D*n_2D;i_bin++)
        RR_2D[i_jack][i_bin]=0;
    }
  }
  for(i_jack=0;i_jack<=n_jack_total;i_jack++){
    for(i_bin=0;i_bin<n_1D;i_bin++){
      DD_1D[i_jack][i_bin] =0;
      DD_l1D[i_jack][i_bin]=0;
      DR_1D[i_jack][i_bin] =0;
      DR_l1D[i_jack][i_bin]=0;
    }
    for(i_bin=0;i_bin<n_2D*n_2D;i_bin++){
      DD_2D[i_jack][i_bin]=0;
      DR_2D[i_jack][i_bin]=0;
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);
  merge_sort(x_data_local,  (size_t)n_data_local,  &x_data_local_index,  SID_REAL,SORT_COMPUTE_INDEX,FALSE);
  merge_sort(x_random_local,(size_t)n_random_local,&x_random_local_index,SID_REAL,SORT_COMPUTE_INDEX,FALSE);
  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
    if(SID.n_proc>1)
      SID_log("Processing rank %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_rank+1,SID.n_proc);
    if(i_rank==0){
      n_data_rank        =n_data_local;
      x_data_rank        =x_data_local;
      y_data_rank        =y_data_local;
      z_data_rank        =z_data_local;
      n_random_rank      =n_random_local;
      x_random_rank      =x_random_local;
      y_random_rank      =y_random_local;
      z_random_rank      =z_random_local;
      data_zone_rank     =data_zone_local;
      random_zone_rank   =random_zone_local;
      x_data_rank_index  =x_data_local_index;
      x_random_rank_index=x_random_local_index;
    }
    // Perform rank exchange
    else{
      if(i_rank==1){
        SID_log("Communicating array sizes...",SID_LOG_OPEN|SID_LOG_CHECKPOINT);
        SID_Allreduce(&n_data_local,  &n_data_allocate,  1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
        SID_log("n_data_max  =%lld",SID_LOG_COMMENT,n_data_allocate);
        SID_Allreduce(&n_random_local,&n_random_allocate,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
        SID_log("n_random_max=%lld",SID_LOG_COMMENT,n_random_allocate);
        SID_log("Done.",SID_LOG_CLOSE);
        SID_log("Allocating arrays...",SID_LOG_OPEN|SID_LOG_CHECKPOINT);
        x_data_rank     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_data_allocate);
        y_data_rank     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_data_allocate);
        z_data_rank     =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_data_allocate);
        x_random_rank   =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_random_allocate);
        y_random_rank   =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_random_allocate);
        z_random_rank   =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_random_allocate);
        data_zone_rank  =(int  *)SID_malloc(sizeof(int)*n_data_allocate);
        random_zone_rank=(int  *)SID_malloc(sizeof(int)*n_random_allocate);        
        SID_log("Done.",SID_LOG_CLOSE);
      }
      SID_log("Performing exchange...",SID_LOG_OPEN|SID_LOG_TIMER);
      exchange_ring_buffer(x_data_local,     sizeof(GBPREAL),n_data_local,  x_data_rank,     &n_data_rank,  i_rank);
      exchange_ring_buffer(y_data_local,     sizeof(GBPREAL),n_data_local,  y_data_rank,     &n_data_rank,  i_rank);
      exchange_ring_buffer(z_data_local,     sizeof(GBPREAL),n_data_local,  z_data_rank,     &n_data_rank,  i_rank);
      exchange_ring_buffer(data_zone_local,  sizeof(GBPREAL),n_data_local,  data_zone_rank,  &n_data_rank,  i_rank);
      exchange_ring_buffer(x_random_local,   sizeof(GBPREAL),n_random_local,x_random_rank,   &n_random_rank,i_rank);
      exchange_ring_buffer(y_random_local,   sizeof(GBPREAL),n_random_local,y_random_rank,   &n_random_rank,i_rank);
      exchange_ring_buffer(z_random_local,   sizeof(GBPREAL),n_random_local,z_random_rank,   &n_random_rank,i_rank);
      exchange_ring_buffer(random_zone_local,sizeof(GBPREAL),n_random_local,random_zone_rank,&n_random_rank,i_rank);
      SID_log("Done.",SID_LOG_CLOSE);
      merge_sort(x_data_rank,  (size_t)n_data_rank,  &x_data_rank_index,  SID_REAL,SORT_COMPUTE_INDEX,FALSE);
      merge_sort(x_random_rank,(size_t)n_random_rank,&x_random_rank_index,SID_REAL,SORT_COMPUTE_INDEX,FALSE);
    }

    // Compute data-data pairs
    SID_log("Processing data-data pairs...",SID_LOG_OPEN|SID_LOG_TIMER|SID_LOG_CHECKPOINT);
    j_data_lo_1=0;             // Limit for the contiguous portion of the periodic box near to the point
    j_data_hi_1=0;             // Limit for the contiguous portion of the periodic box near to the point
    j_data_lo_2=n_data_rank;   // Limit for the (possible) portion of the periodic box at the right edge
    j_data_hi_2=n_data_rank-1; // Limit for the (possible) portion of the periodic box at the right edge
    j_data_lo_1=0;
    j_data_hi_1=n_data_rank-1;
    j_data_lo_2=n_data_rank;
    j_data_hi_2=n_data_rank;
    while((box_size-x_data_rank[x_data_rank_index[j_data_lo_2-1]])<r_max && j_data_lo_2>0){
      j_data_lo_2--;
      if(j_data_lo_2==0) break;
    }
    for(i_data=0;i_data<n_data_local;i_data++){
      x_i=x_data_local[x_data_local_index[i_data]];
      y_i=y_data_local[x_data_local_index[i_data]];
      z_i=z_data_local[x_data_local_index[i_data]];
      // Covers the contiguous range of the periodic box
      while((x_i-x_data_rank[x_data_rank_index[j_data_lo_1]])>r_max && j_data_lo_1<(n_data_rank-2)) j_data_lo_1++;
      while((x_data_rank[x_data_rank_index[j_data_hi_1]]-x_i)<r_max && j_data_hi_1<(n_data_rank-2)) j_data_hi_1++;
      for(j_data=j_data_lo_1;j_data<j_data_hi_1;j_data++){
        if(i_rank!=0 || j_data!=i_data){
          x_j=x_data_rank[x_data_rank_index[j_data]];
          y_j=y_data_rank[x_data_rank_index[j_data]];
          z_j=z_data_rank[x_data_rank_index[j_data]];
          dx =d_periodic(x_i-x_j,box_size);
          dy =d_periodic(y_i-y_j,box_size);
          dz =d_periodic(z_i-z_j,box_size);
          // ... 1D case ...
          sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
          if(sep_1D<r_max_1D){
            bin_l1D=(int)((take_log10(sep_1D)-r_min_l1D)/dr_l1D);
            bin_1D =(int)((sep_1D-r_min_1D)/dr_1D);
            if(bin_1D>=0 && bin_1D<n_1D){
              DD_1D[0][bin_1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(data_zone_local[x_data_local_index[i_data]]!=i_jack && data_zone_rank[x_data_rank_index[j_data]]!=i_jack)
                  DD_1D[i_jack][bin_1D]++;
              }
            }
            if(bin_l1D>=0 && bin_l1D<n_1D){
              DD_l1D[0][bin_l1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(data_zone_local[x_data_local_index[i_data]]!=i_jack && data_zone_rank[x_data_rank_index[j_data]]!=i_jack)
                  DD_l1D[i_jack][bin_l1D]++;
              }
            }
          }
          // ... 2D case ...
          sep_2D_x=sqrt(dx*dx+dy*dy);
          bin_2D_x=(int)((sep_2D_x-r_min_2D)/dr_2D);
          if(bin_2D_x>=0 && bin_2D_x<n_2D){
            sep_2D_y=sqrt(dz*dz);
            bin_2D_y=(int)((sep_2D_y-r_min_2D)/dr_2D);
            if(bin_2D_y>=0 && bin_2D_y<n_2D){
              bin_2D  =bin_2D_y*n_2D+bin_2D_x;
              DD_2D[0][bin_2D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(data_zone_local[x_data_local_index[i_data]]!=i_jack && data_zone_rank[x_data_rank_index[j_data]]!=i_jack)
                  DD_2D[i_jack][bin_2D]++;
              }
            }
          }
        }
      }

      // Covers the (possible) wrap-around end of the periodic box
      if(j_data_lo_2<n_data_rank){
        while((x_i+box_size-x_data_rank[x_data_rank_index[j_data_lo_2]])>r_max){
          j_data_lo_2++;
          if(j_data_lo_2>=n_data_rank)
            break;
        }
      }
      for(j_data=j_data_lo_2;j_data<j_data_hi_2;j_data++){
        if(i_rank!=0 || j_data!=i_data){
          x_j=x_data_rank[x_data_rank_index[j_data]];
          y_j=y_data_rank[x_data_rank_index[j_data]];
          z_j=z_data_rank[x_data_rank_index[j_data]];
          dx =d_periodic(x_i-x_j,box_size);
          dy =d_periodic(y_i-y_j,box_size);
          dz =d_periodic(z_i-z_j,box_size);
          // ... 1D case ...
          sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
          if(sep_1D<r_max_1D){
            bin_l1D=(int)((take_log10(sep_1D)-r_min_l1D)/dr_l1D);
            bin_1D =(int)((sep_1D-r_min_1D)/dr_1D);
            if(bin_1D>=0 && bin_1D<n_1D){
              DD_1D[0][bin_1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(data_zone_local[x_data_local_index[i_data]]!=i_jack && data_zone_rank[x_data_rank_index[j_data]]!=i_jack)
                  DD_1D[i_jack][bin_1D]++;
              }
            }
            if(bin_l1D>=0 && bin_l1D<n_1D){
              DD_l1D[0][bin_l1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(data_zone_local[x_data_local_index[i_data]]!=i_jack && data_zone_rank[x_data_rank_index[j_data]]!=i_jack)
                  DD_l1D[i_jack][bin_l1D]++;
              }
            }
          }
          // ... 2D case ...
          sep_2D_x=sqrt(dx*dx+dy*dy);
          bin_2D_x=(int)((sep_2D_x-r_min_2D)/dr_2D);
          if(bin_2D_x>=0 && bin_2D_x<n_2D){
            sep_2D_y=sqrt(dz*dz);
            bin_2D_y=(int)((sep_2D_y-r_min_2D)/dr_2D);
            if(bin_2D_y>=0 && bin_2D_y<n_2D){
              bin_2D  =bin_2D_y*n_2D+bin_2D_x;
              DD_2D[0][bin_2D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(data_zone_local[x_data_local_index[i_data]]!=i_jack && data_zone_rank[x_data_rank_index[j_data]]!=i_jack)
                  DD_2D[i_jack][bin_2D]++;
              }
            }
          }
        }
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Compute data-random pairs
    SID_log("Processing data-random pairs...",SID_LOG_OPEN|SID_LOG_TIMER|SID_LOG_CHECKPOINT);
    j_random_lo_1=0;             // Limit for the contiguous portion of the periodic box near to the point
    j_random_hi_1=0;             // Limit for the contiguous portion of the periodic box near to the point
    j_random_lo_2=n_random_rank;   // Limit for the (possible) portion of the periodic box at the right edge
    j_random_hi_2=n_random_rank-1; // Limit for the (possible) portion of the periodic box at the right edge
    j_random_lo_1=0;
    j_random_hi_1=n_random_rank-1;
    j_random_lo_2=n_random_rank;
    j_random_hi_2=n_random_rank;
    while((box_size-x_random_rank[x_random_rank_index[j_random_lo_2-1]])<r_max && j_random_lo_2>0){
      j_random_lo_2--;
      if(j_random_lo_2==0) break;
    }
    for(i_data=0;i_data<n_data_local;i_data++){
      x_i=x_data_local[x_data_local_index[i_data]];
      y_i=y_data_local[x_data_local_index[i_data]];
      z_i=z_data_local[x_data_local_index[i_data]];
      // Covers the contiguous range of the periodic box
      while((x_i-x_random_rank[x_random_rank_index[j_random_lo_1]])>r_max && j_random_lo_1<(n_random_rank-2)) j_random_lo_1++;
      while((x_random_rank[x_random_rank_index[j_random_hi_1]]-x_i)<r_max && j_random_hi_1<(n_random_rank-2)) j_random_hi_1++;
      for(j_random=j_random_lo_1;j_random<j_random_hi_1;j_random++){
        x_j=x_random_rank[x_random_rank_index[j_random]];
        y_j=y_random_rank[x_random_rank_index[j_random]];
        z_j=z_random_rank[x_random_rank_index[j_random]];
        dx =d_periodic(x_i-x_j,box_size);
        dy =d_periodic(y_i-y_j,box_size);
        dz =d_periodic(z_i-z_j,box_size);
        // ... 1D case ...
        sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
        if(sep_1D<r_max_1D){
          bin_l1D=(int)((take_log10(sep_1D)-r_min_l1D)/dr_l1D);
          bin_1D =(int)((sep_1D-r_min_1D)/dr_1D);
          if(bin_1D>=0 && bin_1D<n_1D){
            DR_1D[0][bin_1D]++;
            for(i_jack=1;i_jack<=n_jack_total;i_jack++){
              if(data_zone_local[x_data_local_index[i_data]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                DR_1D[i_jack][bin_1D]++;
            }
          }
          if(bin_l1D>=0 && bin_l1D<n_1D){
            DR_l1D[0][bin_l1D]++;
            for(i_jack=1;i_jack<=n_jack_total;i_jack++){
              if(data_zone_local[x_data_local_index[i_data]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                DR_l1D[i_jack][bin_l1D]++;
            }
          }
        }
        // ... 2D case ...
        sep_2D_x=sqrt(dx*dx+dy*dy);
        bin_2D_x=(int)((sep_2D_x-r_min_2D)/dr_2D);
        if(bin_2D_x>=0 && bin_2D_x<n_2D){
          sep_2D_y=sqrt(dz*dz);
          bin_2D_y=(int)((sep_2D_y-r_min_2D)/dr_2D);
          if(bin_2D_y>=0 && bin_2D_y<n_2D){
            bin_2D  =bin_2D_y*n_2D+bin_2D_x;
            DR_2D[0][bin_2D]++;
            for(i_jack=1;i_jack<=n_jack_total;i_jack++){
              if(data_zone_local[x_data_local_index[i_data]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                DR_2D[i_jack][bin_2D]++;
            }
          }
        }
      }

      // Covers the (possible) wrap-around end of the periodic box
      if(j_random_lo_2<n_random_rank){
        while((x_i+box_size-x_random_rank[x_random_rank_index[j_random_lo_2]])>r_max){
          j_random_lo_2++;
          if(j_random_lo_2>=n_random_rank) break;
        }
      }
      for(j_random=j_random_lo_2;j_random<j_random_hi_2;j_random++){
        x_j=x_random_rank[x_random_rank_index[j_random]];
        y_j=y_random_rank[x_random_rank_index[j_random]];
        z_j=z_random_rank[x_random_rank_index[j_random]];
        dx =d_periodic(x_i-x_j,box_size);
        dy =d_periodic(y_i-y_j,box_size);
        dz =d_periodic(z_i-z_j,box_size);
        // ... 1D case ...
        sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
        if(sep_1D<r_max_1D){
          bin_l1D=(int)((take_log10(sep_1D)-r_min_l1D)/dr_l1D);
          bin_1D =(int)((sep_1D-r_min_1D)/dr_1D);
          if(bin_1D>=0 && bin_1D<n_1D){
            DR_1D[0][bin_1D]++;
            for(i_jack=1;i_jack<=n_jack_total;i_jack++){
              if(data_zone_local[x_data_local_index[i_data]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                DR_1D[i_jack][bin_1D]++;
            }
          }
          if(bin_l1D>=0 && bin_l1D<n_1D){
            DR_l1D[0][bin_l1D]++;
            for(i_jack=1;i_jack<=n_jack_total;i_jack++){
              if(data_zone_local[x_data_local_index[i_data]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                DR_l1D[i_jack][bin_l1D]++;
            }
          }
        }
        // ... 2D case ...
        sep_2D_x=sqrt(dx*dx+dy*dy);
        bin_2D_x=(int)((sep_2D_x-r_min_2D)/dr_2D);
        if(bin_2D_x>=0 && bin_2D_x<n_2D){
          sep_2D_y=sqrt(dz*dz);
          bin_2D_y=(int)((sep_2D_y-r_min_2D)/dr_2D);
          if(bin_2D_y>=0 && bin_2D_y<n_2D){
            bin_2D  =bin_2D_y*n_2D+bin_2D_x;
            DR_2D[0][bin_2D]++;
            for(i_jack=1;i_jack<=n_jack_total;i_jack++){
              if(data_zone_local[x_data_local_index[i_data]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                DR_2D[i_jack][bin_2D]++;
            }
          }
        }
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Compute data-random pairs (only done for the first call)
    if((*flag_compute_RR)){
      SID_log("Processing random-random pairs...",SID_LOG_OPEN|SID_LOG_TIMER|SID_LOG_CHECKPOINT);
      j_random_lo_1=0;             // Limit for the contiguous portion of the periodic box near to the point
      j_random_hi_1=0;             // Limit for the contiguous portion of the periodic box near to the point
      j_random_lo_2=n_random_rank;   // Limit for the (possible) portion of the periodic box at the right edge
      j_random_hi_2=n_random_rank-1; // Limit for the (possible) portion of the periodic box at the right edge
      j_random_lo_1=0;
      j_random_hi_1=n_random_rank-1;
      j_random_lo_2=n_random_rank;
      j_random_hi_2=n_random_rank;
      while((box_size-x_random_rank[x_random_rank_index[j_random_lo_2-1]])<r_max && j_random_lo_2>0){
        j_random_lo_2--;
        if(j_random_lo_2==0) break;
      }
      for(i_random=0;i_random<n_random_local;i_random++){
        x_i=x_random_local[x_random_local_index[i_random]];
        y_i=y_random_local[x_random_local_index[i_random]];
        z_i=z_random_local[x_random_local_index[i_random]];
        // Covers the contiguous range of the periodic box
        while((x_i-x_random_rank[x_random_rank_index[j_random_lo_1]])>r_max && j_random_lo_1<(n_random_rank-2)) j_random_lo_1++;
        while((x_random_rank[x_random_rank_index[j_random_hi_1]]-x_i)<r_max && j_random_hi_1<(n_random_rank-2)) j_random_hi_1++;
        for(j_random=j_random_lo_1;j_random<j_random_hi_1;j_random++){
          x_j=x_random_rank[x_random_rank_index[j_random]];
          y_j=y_random_rank[x_random_rank_index[j_random]];
          z_j=z_random_rank[x_random_rank_index[j_random]];
          dx =d_periodic(x_i-x_j,box_size);
          dy =d_periodic(y_i-y_j,box_size);
          dz =d_periodic(z_i-z_j,box_size);
          // ... 1D case ...
          sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
          if(sep_1D<r_max_1D){
            bin_l1D=(int)((take_log10(sep_1D)-r_min_l1D)/dr_l1D);
            bin_1D =(int)((sep_1D-r_min_1D)/dr_1D);
            if(bin_1D>=0 && bin_1D<n_1D){
              RR_1D[0][bin_1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(random_zone_local[x_random_local_index[i_random]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                  RR_1D[i_jack][bin_1D]++;
              }
            }
            if(bin_l1D>=0 && bin_l1D<n_1D){
              RR_l1D[0][bin_l1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(random_zone_local[x_random_local_index[i_random]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                  RR_l1D[i_jack][bin_l1D]++;
              }
            }
          }
          // ... 2D case ...
          sep_2D_x=sqrt(dx*dx+dy*dy);
          bin_2D_x=(int)((sep_2D_x-r_min_2D)/dr_2D);
          if(bin_2D_x>=0 && bin_2D_x<n_2D){
            sep_2D_y=sqrt(dz*dz);
            bin_2D_y=(int)((sep_2D_y-r_min_2D)/dr_2D);
            if(bin_2D_y>=0 && bin_2D_y<n_2D){
              bin_2D  =bin_2D_y*n_2D+bin_2D_x;
              RR_2D[0][bin_2D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(random_zone_local[x_random_local_index[i_random]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                  RR_2D[i_jack][bin_2D]++;
              }
            }
          }
        }

        // Covers the (possible) wrap-around end of the periodic box
        if(j_random_lo_2<n_random_rank){
          while((x_i+box_size-x_random_rank[x_random_rank_index[j_random_lo_2]])>r_max){
            j_random_lo_2++;
            if(j_random_lo_2>=n_random_rank) break;
          }
        }
        for(j_random=j_random_lo_2;j_random<j_random_hi_2;j_random++){
          x_j=x_random_rank[x_random_rank_index[j_random]];
          y_j=y_random_rank[x_random_rank_index[j_random]];
          z_j=z_random_rank[x_random_rank_index[j_random]];
          dx =d_periodic(x_i-x_j,box_size);
          dy =d_periodic(y_i-y_j,box_size);
          dz =d_periodic(z_i-z_j,box_size);
          // ... 1D case ...
          sep_1D=sqrt(dx*dx+dy*dy+dz*dz);
          if(sep_1D<r_max_1D){
            bin_l1D=(int)((take_log10(sep_1D)-r_min_l1D)/dr_l1D);
            bin_1D =(int)((sep_1D-r_min_1D)/dr_1D);
            if(bin_1D>=0 && bin_1D<n_1D){
              RR_1D[0][bin_1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(random_zone_local[x_random_local_index[i_random]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                  RR_1D[i_jack][bin_1D]++;
              }
            }
            if(bin_l1D>=0 && bin_l1D<n_1D){
              RR_l1D[0][bin_l1D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(random_zone_local[x_random_local_index[i_random]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                  RR_l1D[i_jack][bin_l1D]++;
              }
            }
          }
          // ... 2D case ...
          sep_2D_x=sqrt(dx*dx+dy*dy);
          bin_2D_x=(int)((sep_2D_x-r_min_2D)/dr_2D);
          if(bin_2D_x>=0 && bin_2D_x<n_2D){
            sep_2D_y=sqrt(dz*dz);
            bin_2D_y=(int)((sep_2D_y-r_min_2D)/dr_2D);
            if(bin_2D_y>=0 && bin_2D_y<n_2D){
              bin_2D  =bin_2D_y*n_2D+bin_2D_x;
              RR_2D[0][bin_2D]++;
              for(i_jack=1;i_jack<=n_jack_total;i_jack++){
                if(random_zone_local[x_random_local_index[i_random]]!=i_jack && random_zone_rank[x_random_rank_index[j_random]]!=i_jack)
                  RR_2D[i_jack][bin_2D]++;
              }
            }
          }
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);
    }
    if(SID.n_proc>1)
      SID_log("Done.",SID_LOG_CLOSE);
    if(i_rank!=0){
      SID_free(SID_FARG x_data_rank_index);
      SID_free(SID_FARG x_random_rank_index);
    }
  } // i_rank

  if(SID.n_proc>1){
    SID_free(SID_FARG x_data_rank);
    SID_free(SID_FARG y_data_rank);
    SID_free(SID_FARG z_data_rank);
    SID_free(SID_FARG data_zone_rank);
    SID_free(SID_FARG x_random_rank);
    SID_free(SID_FARG y_random_rank);
    SID_free(SID_FARG z_random_rank);
    SID_free(SID_FARG random_zone_rank);
    SID_log("Combining results from separate ranks...",SID_LOG_OPEN);
    for(i_jack=0;i_jack<=n_jack_total;i_jack++){
      SID_Allreduce(SID_IN_PLACE,DD_l1D[i_jack],n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,DR_l1D[i_jack],n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD); 
      SID_Allreduce(SID_IN_PLACE,DD_1D[i_jack], n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,DR_1D[i_jack], n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      if((*flag_compute_RR)){
        SID_Allreduce(SID_IN_PLACE,RR_l1D[i_jack],n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
        SID_Allreduce(SID_IN_PLACE,RR_1D[i_jack], n_1D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      }
      SID_Allreduce(SID_IN_PLACE,DD_2D[i_jack],n_2D*n_2D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,DR_2D[i_jack],n_2D*n_2D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
      if((*flag_compute_RR))
        SID_Allreduce(SID_IN_PLACE,RR_2D[i_jack],n_2D*n_2D,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }

  // Compile correlation functions from data and random pair counts
  SID_log("Processing correlation functions...",SID_LOG_OPEN);
  for(i_bin=0;i_bin<n_1D;i_bin++){
    CFUNC_1D[i_bin] =calc_CFUNC_local((double)DD_1D[0][i_bin]/(double)((n_data)*(n_data-1)),
                                      (double)DR_1D[0][i_bin]/(double)((n_data)*(n_random)),
                                      (double)RR_1D[0][i_bin]/(double)((n_random)*(n_random-1)));
    CFUNC_l1D[i_bin]=calc_CFUNC_local((double)DD_l1D[0][i_bin]/(double)((n_data)*(n_data-1)),
                                      (double)DR_l1D[0][i_bin]/(double)((n_data)*(n_random)),
                                      (double)RR_l1D[0][i_bin]/(double)((n_random)*(n_random-1)));
  }
  for(i_bin=0;i_bin<n_2D*n_2D;i_bin++)
    CFUNC_2D[i_bin]=calc_CFUNC_local((double)DD_2D[0][i_bin]/(double)((n_data)*(n_data-1)),
                                     (double)DR_2D[0][i_bin]/(double)((n_data)*(n_random)),
                                     (double)RR_2D[0][i_bin]/(double)((n_random)*(n_random-1)));
  SID_log("Done.",SID_LOG_CLOSE);

  // Generate 1D covariance matrices
  SID_log("Generating covariance matrices...",SID_LOG_OPEN);
  // ... process log-space profile first ...
  for(i_bin=0;i_bin<n_1D;i_bin++){
    for(i_jack=1,bar_i=0.;i_jack<=n_jack_total;i_jack++){
      bar_i+=calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                              (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                              (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)));
    }
    bar_i/=(double)n_jack_total;
    for(j_bin=0;j_bin<n_1D;j_bin++){
      for(i_jack=1,bar_j=0.,bar_ij=0.;i_jack<=n_jack_total;i_jack++){
        bar_j +=calc_CFUNC_local((double)DD_l1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_l1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_l1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)));
        bar_ij+=calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)))*
                calc_CFUNC_local((double)DD_l1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_l1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_l1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)));
      }
      bar_j /=(double)n_jack_total;
      bar_ij/=(double)n_jack_total;
      COVMTX_l1D[i_bin*n_1D+j_bin]=(double)(n_jack_total-1)*(bar_ij-bar_i*bar_j);
    }
  }
  // ... then process the linear-space profile ...
  for(i_bin=0;i_bin<n_1D;i_bin++){
    for(i_jack=1,bar_i=0.;i_jack<=n_jack_total;i_jack++){
      bar_i+=calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                              (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                              (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)));
    }
    bar_i/=(double)n_jack_total;
    for(j_bin=0;j_bin<n_1D;j_bin++){
      for(i_jack=1,bar_j=0.,bar_ij=0.;i_jack<=n_jack_total;i_jack++){
        bar_j +=calc_CFUNC_local((double)DD_1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)));
        bar_ij+=calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)))*
                calc_CFUNC_local((double)DD_1D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_1D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_1D[i_jack][j_bin]/(double)((n_random)*(n_random-1)));
      }
      bar_j /=(double)n_jack_total;
      bar_ij/=(double)n_jack_total;
      COVMTX_1D[i_bin*n_1D+j_bin]=(double)(n_jack_total-1)*(bar_ij-bar_i*bar_j);
    }
  }
  // ... and finally, the 2D case ...
  for(i_bin=0;i_bin<n_2D_total;i_bin++){
    for(i_jack=1,bar_i=0.;i_jack<=n_jack_total;i_jack++){
      bar_i+=calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                              (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                              (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)));
    }
    bar_i/=(double)n_jack_total;
    for(j_bin=0;j_bin<n_2D_total;j_bin++){
      for(i_jack=1,bar_j=0.,bar_ij=0.;i_jack<=n_jack_total;i_jack++){
        bar_j +=calc_CFUNC_local((double)DD_2D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_2D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_2D[i_jack][j_bin]/(double)((n_random)*(n_random-1)));
        bar_ij+=calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)))*
                calc_CFUNC_local((double)DD_2D[i_jack][j_bin]/(double)((n_data)*(n_data-1)),
                                 (double)DR_2D[i_jack][j_bin]/(double)((n_data)*(n_random)),
                                 (double)RR_2D[i_jack][j_bin]/(double)((n_random)*(n_random-1)));
      }
      bar_j /=(double)n_jack_total;
      bar_ij/=(double)n_jack_total;
      COVMTX_2D[i_bin*n_2D_total+j_bin]=(double)(n_jack_total-1)*(bar_ij-bar_i*bar_j);
    }
  }
  SID_log("Done.",SID_LOG_CLOSE);


  // Compile uncertainties from jacknife-resampling data and random pair counts
  SID_log("Processing jack knifes...",SID_LOG_OPEN);
  // ... 1D case ...
  for(i_bin=0;i_bin<n_1D;i_bin++){
    for(i_jack=1,mean=0.;i_jack<=n_jack_total;i_jack++)
      mean+=calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                             (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                             (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)));    
    mean/=(double)n_jack_total;
    for(i_jack=1,std_dev=0.;i_jack<=n_jack_total;i_jack++)
      std_dev+=pow(mean-calc_CFUNC_local((double)DD_l1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                         (double)DR_l1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                         (double)RR_l1D[i_jack][i_bin]/(double)((n_random)*(n_random-1))),2.);
    std_dev=sqrt(std_dev/(double)n_jack_total);
    dCFUNC_l1D[i_bin]=std_dev*sqrt((double)(n_jack_total-1));
    for(i_jack=1,mean=0.;i_jack<=n_jack_total;i_jack++)
      mean+=calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                             (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                             (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1)));    
    mean/=(double)n_jack_total;
    for(i_jack=1,std_dev=0.;i_jack<=n_jack_total;i_jack++)
      std_dev+=pow(mean-calc_CFUNC_local((double)DD_1D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                         (double)DR_1D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                         (double)RR_1D[i_jack][i_bin]/(double)((n_random)*(n_random-1))),2.);
    std_dev=sqrt(std_dev/(double)n_jack_total);
    dCFUNC_1D[i_bin]=std_dev*sqrt((double)(n_jack_total-1));
  }

  // ... 2D case ...
  for(i_bin=0;i_bin<n_2D*n_2D;i_bin++){
    for(i_jack=1,mean=0.;i_jack<=n_jack_total;i_jack++)
      mean+=calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                             (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                             (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1)));    
    mean/=(double)n_jack_total;
    for(i_jack=1,std_dev=0.;i_jack<=n_jack_total;i_jack++)
      std_dev+=pow(mean-calc_CFUNC_local((double)DD_2D[i_jack][i_bin]/(double)((n_data)*(n_data-1)),
                                         (double)DR_2D[i_jack][i_bin]/(double)((n_data)*(n_random)),
                                         (double)RR_2D[i_jack][i_bin]/(double)((n_random)*(n_random-1))),2.);
    std_dev=sqrt(std_dev/(double)n_jack_total);
    dCFUNC_2D[i_bin]=std_dev*sqrt((double)(n_jack_total-1));
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  SID_free(SID_FARG data_zone_local);
  SID_free(SID_FARG random_zone_local);
  SID_free(SID_FARG x_data_local_index);
  SID_free(SID_FARG x_random_local_index);
  if(flag_alloc_x)
    SID_free(SID_FARG x_data_local);
  if(flag_alloc_y)
    SID_free(SID_FARG y_data_local);
  if(flag_alloc_z)
    SID_free(SID_FARG z_data_local);
  SID_log("Done.",SID_LOG_CLOSE);
  (*flag_compute_RR)=FALSE;  
  SID_log("Done.",SID_LOG_CLOSE);
}

