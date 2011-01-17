#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

void cell2index(int cell[3],int n_grid[3],size_t *index){
  int i_x,i_y,i_z;
  i_x=cell[0];
  i_y=cell[1];
  i_z=cell[2];
  if(i_x<0)
    i_x+=n_grid[0];
  else if(i_x>=n_grid[0])
    i_x-=n_grid[0];
  if(i_y<0)
    i_y+=n_grid[1];
  else if(i_y>=n_grid[1])
    i_y-=n_grid[1];
  if(i_z<0)
    i_z+=n_grid[2];
  else if(i_z>=n_grid[2])
    i_z-=n_grid[2];
  (*index)=(i_x*n_grid[1]+i_y)*n_grid[2]+i_z;
}

void  create_local_buffer(int     n_grid[3],
                          int     x_i_cell,
                          int     y_i_cell,
                          int     z_i_cell,
                          int     i_cell,
                          int     buffer_size,
                          int     flag_periodic_box,
                          int     flag_use_velocities,
                          int     flag_use_mass,
                          int     flag_multimass,
                          size_t *cell_index_local,
                          size_t *cell_index_sort,
                          size_t  n_particles_local,
                          REAL   *x_particles_local,
                          REAL   *y_particles_local,
                          REAL   *z_particles_local,
                          REAL   *vx_particles_local,
                          REAL   *vy_particles_local,
                          REAL   *vz_particles_local,
                          double *m_particles_local,
                          REAL   *x_particles_buffer,
                          REAL   *y_particles_buffer,
                          REAL   *z_particles_buffer,
                          REAL   *vx_particles_buffer,
                          REAL   *vy_particles_buffer,
                          REAL   *vz_particles_buffer,
                          double *m_particles_buffer,
                          size_t *n_particles_buffer,
                          int     n_buffer_max);
void  create_local_buffer(int     n_grid[3],
                          int     x_i_cell,
                          int     y_i_cell,
                          int     z_i_cell,
                          int     i_cell,
                          int     buffer_size,
                          int     flag_periodic_box,
                          int     flag_use_velocities,
                          int     flag_use_mass,
                          int     flag_multimass,
                          size_t *cell_index_local,
                          size_t *cell_index_sort,
                          size_t  n_particles_local,
                          REAL   *x_particles_local,
                          REAL   *y_particles_local,
                          REAL   *z_particles_local,
                          REAL   *vx_particles_local,
                          REAL   *vy_particles_local,
                          REAL   *vz_particles_local,
                          double *m_particles_local,
                          REAL   *x_particles_buffer,
                          REAL   *y_particles_buffer,
                          REAL   *z_particles_buffer,
                          REAL   *vx_particles_buffer,
                          REAL   *vy_particles_buffer,
                          REAL   *vz_particles_buffer,
                          double *m_particles_buffer,
                          size_t *n_particles_buffer,
                          int     n_buffer_max){
  size_t i_particle;
  size_t j_particle;
  size_t start_particle_buffer;
  size_t n_buffer;
  size_t buffer_cell_index;
  int    buffer_cell_x;
  int    buffer_cell_y;
  int    buffer_cell_z;
  int    buffer_cell[3];

  // Sweep buffer cells in x
  (*n_particles_buffer)=0;
  n_buffer=0;
  for(buffer_cell_x=x_i_cell-buffer_size;buffer_cell_x<=x_i_cell+buffer_size;buffer_cell_x++){
    buffer_cell[0]=buffer_cell_x;
    // Deal with periodic/non-periodic BCs in x
    if(flag_periodic_box && buffer_cell[0]<0)          buffer_cell[0]+=n_grid[0];
    if(flag_periodic_box && buffer_cell[0]>=n_grid[0]) buffer_cell[0]-=n_grid[0];
    if(buffer_cell[0]>=0 && buffer_cell[0]<n_grid[0]){
      // Sweep buffer cells in y
      for(buffer_cell_y=y_i_cell-buffer_size;buffer_cell_y<=y_i_cell+buffer_size;buffer_cell_y++){
        buffer_cell[1]=buffer_cell_y;
        // Deal with periodic/non-periodic BCs in y
        if(flag_periodic_box && buffer_cell[1]<0)          buffer_cell[1]+=n_grid[1];
        if(flag_periodic_box && buffer_cell[1]>=n_grid[1]) buffer_cell[1]-=n_grid[1];
        if(buffer_cell[1]>=0 && buffer_cell[1]<n_grid[1]){
          // Sweep buffer cells in z
          for(buffer_cell_z=z_i_cell-buffer_size;buffer_cell_z<=z_i_cell+buffer_size;buffer_cell_z++){
            buffer_cell[2]=buffer_cell_z;
            // Deal with periodic/non-periodic BCs in z
            if(flag_periodic_box && buffer_cell[2]<0)          buffer_cell[2]+=n_grid[2];
            if(flag_periodic_box && buffer_cell[2]>=n_grid[2]) buffer_cell[2]-=n_grid[2];
            if(buffer_cell[2]>=0 && buffer_cell[2]<n_grid[2]){
              // Only enter cells that are on the outer shell of the buffer
              if(buffer_cell_x==x_i_cell-buffer_size ||
                 buffer_cell_x==x_i_cell+buffer_size ||
                 buffer_cell_y==y_i_cell-buffer_size ||
                 buffer_cell_y==y_i_cell+buffer_size ||
                 buffer_cell_z==z_i_cell-buffer_size ||
                 buffer_cell_z==z_i_cell+buffer_size){
                // Count particles in this buffer cell;
                //   j_particle runs over all the (sorted) particles in the current buffer cell
                cell2index(buffer_cell,n_grid,&buffer_cell_index);
                start_particle_buffer=find_index(cell_index_local,buffer_cell_index,n_particles_local,cell_index_sort);
                j_particle           =start_particle_buffer;
                while(cell_index_local[cell_index_sort[j_particle]]==buffer_cell_index){
                  if(n_buffer>=(size_t)n_buffer_max)
                    SID_trap_error("Buffer size too small in create_local_buffer",ERROR_LOGIC);
                  x_particles_buffer[n_buffer]=x_particles_local[cell_index_sort[j_particle]];
                  y_particles_buffer[n_buffer]=y_particles_local[cell_index_sort[j_particle]];
                  z_particles_buffer[n_buffer]=z_particles_local[cell_index_sort[j_particle]];
                  if(flag_use_velocities){
                    vx_particles_buffer[n_buffer]=vx_particles_local[cell_index_sort[j_particle]];
                    vy_particles_buffer[n_buffer]=vy_particles_local[cell_index_sort[j_particle]];
                    vz_particles_buffer[n_buffer]=vz_particles_local[cell_index_sort[j_particle]];                    
                  }
                  if(flag_use_mass && flag_multimass)
                    m_particles_buffer[n_buffer]=m_particles_local[cell_index_sort[j_particle]];                                        
                  n_buffer++;
                  j_particle++;
                  if(j_particle==n_particles_local) break;
                }
              } // if cell is on buffer shell
            } // force periodic in z
          } // sweep cells in z
        } // force periodic in y
      } // sweep cells in y
    } // force periodic in x
  } // sweep cells in x
  (*n_particles_buffer)=n_buffer;
}

void  exchange_buffer(int     i_rank,
                      int     x_i_cell_requested,
                      int     y_i_cell_requested,
                      int     z_i_cell_requested,
                      int     i_cell_requested,
                      int     buffer_size_requested,
                      int     n_grid[3],
                      int     flag_periodic_box,
                      int     flag_use_velocities,
                      int     flag_use_mass,
                      int     flag_multimass,
                      size_t *cell_index_local,
                      size_t *cell_index_sort,
                      size_t  n_particles_local,
                      REAL   *x_particles_local,
                      REAL   *y_particles_local,
                      REAL   *z_particles_local,
                      REAL   *vx_particles_local,
                      REAL   *vy_particles_local,
                      REAL   *vz_particles_local,
                      double *m_particles_local,
                      REAL   *x_particles_buffer,
                      REAL   *y_particles_buffer,
                      REAL   *z_particles_buffer,
                      REAL   *vx_particles_buffer,
                      REAL   *vy_particles_buffer,
                      REAL   *vz_particles_buffer,
                      double *m_particles_buffer,
                      size_t *n_particles_buffer,
                      int     n_buffer_max);
void  exchange_buffer(int     i_rank,
                      int     x_i_cell_requested,
                      int     y_i_cell_requested,
                      int     z_i_cell_requested,
                      int     i_cell_requested,
                      int     buffer_size_requested,
                      int     n_grid[3],
                      int     flag_periodic_box,
                      int     flag_use_velocities,
                      int     flag_use_mass,
                      int     flag_multimass,
                      size_t *cell_index_local,
                      size_t *cell_index_sort,
                      size_t  n_particles_local,
                      REAL   *x_particles_local,
                      REAL   *y_particles_local,
                      REAL   *z_particles_local,
                      REAL   *vx_particles_local,
                      REAL   *vy_particles_local,
                      REAL   *vz_particles_local,
                      double *m_particles_local,
                      REAL   *x_particles_buffer,
                      REAL   *y_particles_buffer,
                      REAL   *z_particles_buffer,
                      REAL   *vx_particles_buffer,
                      REAL   *vy_particles_buffer,
                      REAL   *vz_particles_buffer,
                      double *m_particles_buffer,
                      size_t *n_particles_buffer,
                      int     n_buffer_max){
  size_t  i_particle;
  size_t  n_particles_send;
  void   *send_buffer;
  int     x_i_cell;
  int     y_i_cell;
  int     z_i_cell;
  int     i_cell;
  int     buffer_size;
  int     rank_to,rank_from;

  // Decide which ranks are swapping with which
  rank_to  =(SID.My_rank+i_rank);
  if(rank_to>=SID.n_proc)
    rank_to-=SID.n_proc;
  rank_from=(SID.My_rank-i_rank);
  if(rank_from<0)
    rank_from+=SID.n_proc;

#if USE_MPI
  if(rank_to!=rank_from){
    MPI_Sendrecv(&x_i_cell_requested,
                 1,
                 MPI_INTEGER,
                 rank_to,
                 121,
                 &x_i_cell,
                 1,
                 MPI_INTEGER,
                 rank_from,
                 121,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(&y_i_cell_requested,
                 1,
                 MPI_INTEGER,
                 rank_to,
                 122,
                 &y_i_cell,
                 1,
                 MPI_INTEGER,
                 rank_from,
                 122,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(&z_i_cell_requested,
                 1,
                 MPI_INTEGER,
                 rank_to,
                 123,
                 &z_i_cell,
                 1,
                 MPI_INTEGER,
                 rank_from,
                 123,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(&i_cell_requested,
                 1,
                 MPI_INTEGER,
                 rank_to,
                 124,
                 &i_cell,
                 1,
                 MPI_INTEGER,
                 rank_from,
                 124,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(&buffer_size_requested,
                 1,
                 MPI_INTEGER,
                 rank_to,
                 125,
                 &buffer_size,
                 1,
                 MPI_INTEGER,
                 rank_from,
                 125,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
  }
  else{
    x_i_cell   =x_i_cell_requested;
    y_i_cell   =y_i_cell_requested;
    z_i_cell   =z_i_cell_requested;
    i_cell     =i_cell_requested;
    buffer_size=buffer_size_requested;    
  }
#else
  x_i_cell   =x_i_cell_requested;
  y_i_cell   =y_i_cell_requested;
  z_i_cell   =z_i_cell_requested;
  i_cell     =i_cell_requested;
  buffer_size=buffer_size_requested;
#endif

  create_local_buffer(n_grid,
                      x_i_cell,
                      y_i_cell,
                      z_i_cell,
                      i_cell,
                      buffer_size,
                      flag_periodic_box,
                      flag_use_velocities,
                      flag_use_mass,
                      flag_multimass,
                      cell_index_local,
                      cell_index_sort,
                      n_particles_local,
                      x_particles_local,
                      y_particles_local,
                      z_particles_local,
                      vx_particles_local,
                      vy_particles_local,
                      vz_particles_local,
                      m_particles_local,
                      x_particles_buffer,
                      y_particles_buffer,
                      z_particles_buffer,
                      vx_particles_buffer,
                      vy_particles_buffer,
                      vz_particles_buffer,
                      m_particles_buffer,
                      &n_particles_send,
                      n_buffer_max);

#ifdef USE_MPI
  // Exchange this buffer cell's particles with other ranks
  if(rank_to!=rank_from){

    // Create send buffer
    if(!(flag_use_mass && flag_multimass))
      send_buffer=SID_malloc(n_particles_send*sizeof(REAL));
    else
      send_buffer=SID_malloc(n_particles_send*sizeof(double));

    // Exchange number of particles to exchange ...
    MPI_Sendrecv(&n_particles_send,
                 1,
                 MPI_SIZE_T,
                 rank_to,
                 126,
                 n_particles_buffer,
                 1,
                 MPI_SIZE_T,
                 rank_from,
                 126,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    // ... positions ...
    for(i_particle=0;i_particle<n_particles_send;i_particle++)
      ((REAL *)send_buffer)[i_particle]=x_particles_buffer[i_particle];
    MPI_Sendrecv(send_buffer,
                 (int)n_particles_send,
                 MPI_REAL,
                 rank_to,
                 127,
                 x_particles_buffer,
                 (int)(*n_particles_buffer),
                 MPI_REAL,
                 rank_from,
                 127,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    for(i_particle=0;i_particle<n_particles_send;i_particle++)
      ((REAL *)send_buffer)[i_particle]=y_particles_buffer[i_particle];
    MPI_Sendrecv(send_buffer,
                 (int)n_particles_send,
                 MPI_REAL,
                 rank_to,
                 128,
                 y_particles_buffer,
                 (int)(*n_particles_buffer),
                 MPI_REAL,
                 rank_from,
                 128,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    for(i_particle=0;i_particle<n_particles_send;i_particle++)
      ((REAL *)send_buffer)[i_particle]=z_particles_buffer[i_particle];
    MPI_Sendrecv(send_buffer,
                 (int)n_particles_send,
                 MPI_REAL,
                 rank_to,
                 129,
                 z_particles_buffer,
                 (int)(*n_particles_buffer),
                 MPI_REAL,
                 rank_from,
                 129,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    // ... velocities ...
    if(flag_use_velocities){
      for(i_particle=0;i_particle<n_particles_send;i_particle++)
        ((REAL *)send_buffer)[i_particle]=vx_particles_buffer[i_particle];
      MPI_Sendrecv(send_buffer,
                   (int)n_particles_send,
                   MPI_REAL,
                   rank_to,
                   130,
                   vx_particles_buffer,
                   (int)(*n_particles_buffer),
                   MPI_REAL,
                   rank_from,
                   130,
                   MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
      for(i_particle=0;i_particle<n_particles_send;i_particle++)
        ((REAL *)send_buffer)[i_particle]=vy_particles_buffer[i_particle];
      MPI_Sendrecv(send_buffer,
                   (int)n_particles_send,
                   MPI_REAL,
                   rank_to,
                   131,
                   vy_particles_buffer,
                   (int)(*n_particles_buffer),
                   MPI_REAL,
                   rank_from,
                   131,
                   MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
      for(i_particle=0;i_particle<n_particles_send;i_particle++)
        ((REAL *)send_buffer)[i_particle]=vz_particles_buffer[i_particle];
      MPI_Sendrecv(send_buffer,
                   (int)n_particles_send,
                   MPI_REAL,
                   rank_to,
                   132,
                   vz_particles_buffer,
                   (int)(*n_particles_buffer),
                   MPI_REAL,
                   rank_from,
                   132,
                   MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
    }
    if(flag_use_mass && flag_multimass){
      for(i_particle=0;i_particle<n_particles_send;i_particle++)
        ((REAL *)send_buffer)[i_particle]=m_particles_buffer[i_particle];
      MPI_Sendrecv(send_buffer,
                   (int)n_particles_send,
                   MPI_REAL,
                   rank_to,
                   133,
                   m_particles_buffer,
                   (int)(*n_particles_buffer),
                   MPI_REAL,
                   rank_from,
                   133,
                   MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
  
    }
    SID_free((void **)&send_buffer);
  }
  else
    (*n_particles_buffer)=n_particles_send;
#else
  (*n_particles_buffer)=n_particles_send;
#endif
}

void smooth_grid(plist_info *plist,
                 int         grid_size,
                 char       *species_name,
                 int         mode){
  size_t      i_p;
  size_t      i_cell,j_cell,k_cell;
  size_t      n_p_cell;
  size_t      n_p_cell_local;
  size_t      n_p_cell_max;
  size_t      n_p_cell_max_local;
  size_t      n_particles_buffer;
  size_t      n_cells;
  size_t      start_particle;
  size_t      start_particle_buffer;
  int         flag_exchange;
  int         x_i_cell;
  int         y_i_cell;
  int         z_i_cell;
  int         buffer_cell[3];
  int         buffer_cell_x;
  int         buffer_cell_y;
  int         buffer_cell_z;
  size_t      buffer_cell_index;
  size_t      n_p_buffer_cell;
  char        var_name[256];
  size_t      n_particles_local;
  size_t      n_particles;
  size_t     *cell_index_local;
  REAL       *x_particles_local;
  REAL       *y_particles_local;
  REAL       *z_particles_local;
  REAL       *vx_particles_local;
  REAL       *vy_particles_local;
  REAL       *vz_particles_local;
  double     *m_particles_local;
  REAL       *x_particles_buffer;
  REAL       *y_particles_buffer;
  REAL       *z_particles_buffer;
  REAL       *vx_particles_buffer;
  REAL       *vy_particles_buffer;
  REAL       *vz_particles_buffer;
  double     *m_particles_buffer;
  void       *send_buffer;
  double      m_p;
  int         flag_multimass;
  size_t      n_send;
  size_t      send_size;
  size_t      receive_left_size=0;
  size_t      receive_right_size=0;
  REAL       *send_left;
  REAL       *send_right;
  REAL       *receive_left=NULL;
  REAL       *receive_right=NULL;
  double       r_i,r_min;
  double       box_size;
  double       W_i;
  double       L_grid[3];
  int          n_grid[3];
  size_t       cells_per_slice;
  size_t      *cell_index_sort,*i_p_slice_start,*i_p_slice_stop,*n_p_slice,i_particle;
  size_t       n_left;
  size_t       n_left_cell_local;
  size_t       i_slice,N_neighbours;
  REAL        *r_smooth;
  REAL        *rho_smooth;
  double       rho_temp;
  REAL        *sigma_v_smooth;
  int          i_neighbour;
  REAL       **r_neighbours;
  REAL       **v_neighbours;
  double     **m_neighbours;
  short int   *n_neighbours;
  int         *status;
  size_t       j_particle;
  int          k_particle;
  int          l_particle;
  int          i_table;
  int          buffer_size;
  int          buffer_size_max;
  REAL         r_ij;
  int          flag_break_early;
  int          flag_use_mass;
  int          flag_use_velocities;
  int          flag_calc_rho;
  int          flag_calc_sigma_v;
  int          flag_periodic_box;
  int          i_rank;
  int          rank_to;
  int          rank_from;
  int          n_buffer_max;
  int          flag_rank_alive;
  int          n_ranks_left;
  size_t       n_left_total_local;
  size_t       n_left_total;
  int          n_update;
  size_t       n_next_update;
  size_t       n_update_step;
  time_t       time_elapsed,time_start,time_remaining;
  char         time_string[50];
  double      *kernel_radius;
  double      *kernel_table;
  double       f_table;
  double       w_kernel;
  double       inv_r_smooth3;
  

  SID_log("Computing smoothed quantities for %s particles...",SID_LOG_OPEN|SID_LOG_TIMER,species_name);

  // Treat the box as periodic? (defualt is yes)
  if(check_mode_for_flag(mode,SMOOTH_NOT_PERIODIC))
    flag_periodic_box=FALSE;
  else
    flag_periodic_box=TRUE;

  // Generate sph kernel
  set_sph_kernel(plist,SPH_KERNEL_GADGET);
  kernel_radius=(double *)ADaPS_fetch(plist->data,"sph_kernel_radius");
  kernel_table =(double *)ADaPS_fetch(plist->data,"sph_kernel_3d");  

  // Set grid
  box_size=((double *)ADaPS_fetch(plist->data,"box_size"))[0];
  L_grid[0]=box_size;
  L_grid[1]=box_size;
  L_grid[2]=box_size;
  n_grid[0]=grid_size;
  n_grid[1]=grid_size;
  n_grid[2]=grid_size;
  n_cells  =
    (size_t)n_grid[0]*
    (size_t)n_grid[1]*
    (size_t)n_grid[2];
  cells_per_slice =
    (size_t)n_grid[1]*
    (size_t)n_grid[2];

  // Fetch the needed particle information
  n_particles        =((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0]; 
  n_particles_local  =((size_t *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0]; 
  x_particles_local  = (REAL   *)ADaPS_fetch(plist->data,"x_%s",    species_name);
  y_particles_local  = (REAL   *)ADaPS_fetch(plist->data,"y_%s",    species_name);
  z_particles_local  = (REAL   *)ADaPS_fetch(plist->data,"z_%s",    species_name);

  // Assign particles to the grid
  SID_log("Assigning %s particles to grid cells...",SID_LOG_OPEN|SID_LOG_TIMER,species_name);
  assign_particles_to_grid(plist,L_grid,n_grid,species_name);
  cell_index_local=(size_t *)ADaPS_fetch(plist->data,"cell_index_%s",species_name);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Sorting grid cell indices...",SID_LOG_OPEN|SID_LOG_TIMER,species_name);
  sort(cell_index_local,n_particles_local,&cell_index_sort,ADaPS_SIZE_T,SORT_LOCAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  SID_log("Done.",SID_LOG_CLOSE);

  // Determine particle list offsets and sizes for each y-z slice
  i_p_slice_start=(size_t *)SID_malloc(sizeof(size_t)*grid_size);
  i_p_slice_stop =(size_t *)SID_malloc(sizeof(size_t)*grid_size);
  n_p_slice      =(size_t *)SID_malloc(sizeof(size_t)*grid_size);
  for(i_slice=0,i_particle=0,flag_break_early=FALSE;i_slice<grid_size && !flag_break_early;i_slice++){
    i_p_slice_start[i_slice]=i_particle;
    i_p_slice_stop[i_slice] =i_particle;
    n_p_slice[i_slice]      =0;
    while(i_slice==cell_index_local[cell_index_sort[i_particle]]/cells_per_slice){
      i_p_slice_stop[i_slice]=i_particle++;
      n_p_slice[i_slice]++;
      if(i_particle==n_particles_local){
        flag_break_early=TRUE;
        break;
      }
    }
  }

  // Allocate arrays for storing final results and fetch some (optional) needed particle info
  N_neighbours     =64;
  flag_use_mass    =FALSE;
  flag_calc_rho    =FALSE;
  flag_calc_sigma_v=FALSE;
  r_smooth         =(REAL *)SID_malloc(sizeof(REAL)*n_particles_local);
  if(check_mode_for_flag(mode,SMOOTH_DENSITY) || check_mode_for_flag(mode,SMOOTH_DEFAULT)){
    rho_smooth   =(REAL *)SID_malloc(sizeof(REAL)*n_particles_local);
    flag_calc_rho=TRUE;
    flag_use_mass=TRUE;
  }
  if(check_mode_for_flag(mode,SMOOTH_SIGMA_V) || check_mode_for_flag(mode,SMOOTH_DEFAULT)){
    vx_particles_local = (REAL *)ADaPS_fetch(plist->data,"vx_%s",species_name);
    vy_particles_local = (REAL *)ADaPS_fetch(plist->data,"vy_%s",species_name);
    vz_particles_local = (REAL *)ADaPS_fetch(plist->data,"vz_%s",species_name);
    sigma_v_smooth     = (REAL *)SID_malloc(sizeof(REAL)*n_particles_local);
    flag_calc_sigma_v  = TRUE;
    flag_use_velocities= TRUE;
    flag_use_mass      = TRUE;
  }
  if(flag_use_mass){
    if(ADaPS_exist(plist->data,"M_%s",species_name)){
      m_particles_local=(double *)ADaPS_fetch(plist->data,"M_%s",species_name);
      flag_multimass=TRUE;
    }
    else{
      if(ADaPS_exist(plist->data,"mass_array_%s",species_name))
        m_p=((double *)ADaPS_fetch(plist->data,"mass_array_%s",species_name))[0];
      else
        m_p=1.;
      flag_multimass=FALSE;
    }
  }
  
  // Determine how big the neighbour arrays need to be
  for(i_cell=0,i_particle=0,n_p_cell_max_local=0,flag_break_early=FALSE;i_cell<n_cells && !flag_break_early;i_cell++){
    n_p_cell=0;
    while(cell_index_local[cell_index_sort[i_particle]]==i_cell){
      n_p_cell++;
      i_particle++;
      if(i_particle==n_particles_local){
        flag_break_early=TRUE;
        break;
      }
    }
    n_p_cell_max_local=MAX(n_p_cell_max_local,n_p_cell);
  }
#ifdef USE_MPI
  MPI_Allreduce(&n_p_cell_max_local,&n_p_cell_max,1,MPI_SIZE_T,MPI_MAX,MPI_COMM_WORLD);
#else
  n_p_cell_max=n_p_cell_max_local;
#endif

  // Allocate arrays for storing neighbour info
  n_neighbours    =(short int  *)SID_malloc(sizeof(short int)*n_p_cell_max);
  status          =(int        *)SID_malloc(sizeof(int      )*n_p_cell_max);
  r_neighbours    =(REAL      **)SID_malloc(sizeof(REAL    *)*n_p_cell_max);
  if(flag_use_velocities)
    v_neighbours=(REAL **)SID_malloc(sizeof(REAL *)*n_p_cell_max);
  if(flag_use_mass && flag_multimass)
    m_neighbours=(double **)SID_malloc(sizeof(double *)*n_p_cell_max);
  for(k_particle=0;k_particle<n_p_cell_max;k_particle++){
    r_neighbours[k_particle] =(REAL   *)SID_malloc(sizeof(REAL)  *N_neighbours);
    if(flag_use_velocities)
      v_neighbours[k_particle]=(REAL *)SID_malloc(sizeof(REAL)*N_neighbours);
    if(flag_use_mass && flag_multimass)
      m_neighbours[k_particle]=(double *)SID_malloc(sizeof(double)*N_neighbours);
  }

  // Allocate buffer arrays
  n_buffer_max=100*n_p_cell_max;
  x_particles_buffer =(REAL   *)SID_malloc(sizeof(REAL)*n_buffer_max);
  y_particles_buffer =(REAL   *)SID_malloc(sizeof(REAL)*n_buffer_max);
  z_particles_buffer =(REAL   *)SID_malloc(sizeof(REAL)*n_buffer_max);
  if(flag_use_velocities){
    vx_particles_buffer=(REAL   *)SID_malloc(sizeof(REAL)*n_buffer_max);
    vy_particles_buffer=(REAL   *)SID_malloc(sizeof(REAL)*n_buffer_max);
    vz_particles_buffer=(REAL   *)SID_malloc(sizeof(REAL)*n_buffer_max);
  }
  if(flag_use_mass && flag_multimass)
    m_particles_buffer =(double *)SID_malloc(sizeof(double)*n_buffer_max);

  // Set some stuff for the run-time status messages
  n_update     =100;
  n_update_step=(size_t)((double)n_particles/(double)n_update);
  n_next_update=n_particles-n_update_step;
  time(&time_start);

  // Process particles slice-by-slice
  flag_rank_alive   =1;
  n_ranks_left      =SID.n_proc;
  n_left_total_local=n_particles_local;
  n_left_total      =n_particles;
  while(n_ranks_left>0 && n_left_total>0){
    for(x_i_cell=0,i_cell=0,i_slice=0;x_i_cell<n_grid[0] && flag_rank_alive;x_i_cell++,i_slice++){
      // Process slice cell-by-cell (they have a common buffer)
      buffer_size_max=-1; // keep track of largest buffer used for this slice
      for(y_i_cell=0;y_i_cell<n_grid[1];y_i_cell++){
        for(z_i_cell=0;z_i_cell<n_grid[2];z_i_cell++,i_cell++){
          // Count particles in the cell currently being smoothed and initialize arrays
          start_particle=find_index(cell_index_local,i_cell,n_particles_local,cell_index_sort);
          j_particle    =start_particle;
          k_particle    =0; // Runs over the (sorted) particles in the smoothing cell; corresponds to i_particle which runs over all local (sorted) particles
          n_p_cell_local=0;
          if(j_particle<n_particles_local){
            while(cell_index_local[cell_index_sort[j_particle]]==i_cell){
              n_neighbours[k_particle]     =0;
              status[k_particle]           =4*grid_size;
              k_particle++;
              n_p_cell_local++;
              j_particle++;
              if(j_particle==n_particles_local) break;
            }
          }
          // Keep expanding the buffer until all particles are done
          n_left_cell_local=n_p_cell_local;
          buffer_size =-1;
          while(n_left_cell_local>0){
            buffer_size++;
            buffer_size_max=MAX(buffer_size_max,buffer_size);
            // Check buffer particles from each rank in turn
#ifdef USE_MPI
            MPI_Allreduce(&flag_rank_alive,&n_ranks_left,   1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(&n_left_total_local,&n_left_total,1,MPI_SIZE_T, MPI_SUM,MPI_COMM_WORLD);
#else
            n_ranks_left=flag_rank_alive;
            n_left_total=n_left_total_local;
#endif
            for(i_rank=0;i_rank<SID.n_proc;i_rank++){
              exchange_buffer(i_rank,
                              x_i_cell,
                              y_i_cell,
                              z_i_cell,
                              i_cell,
                              buffer_size,
                              n_grid,
                              flag_periodic_box,
                              flag_use_velocities,
                              flag_use_mass,
                              flag_multimass,
                              cell_index_local,
                              cell_index_sort,
                              n_particles_local,
                              x_particles_local,
                              y_particles_local,
                              z_particles_local,
                              vx_particles_local,
                              vy_particles_local,
                              vz_particles_local,
                              m_particles_local,
                              x_particles_buffer,
                              y_particles_buffer,
                              z_particles_buffer,
                              vx_particles_buffer,
                              vy_particles_buffer,
                              vz_particles_buffer,
                              m_particles_buffer,
                              &n_particles_buffer,
                              n_buffer_max);
              // Scan over buffer particles
              for(k_particle=0,i_particle=start_particle;k_particle<n_p_cell_local;i_particle++,k_particle++){
                if(status[k_particle]>0){
                  // Perform scan on buffer particles
                  for(j_particle=0;j_particle<n_particles_buffer;j_particle++){
                    if(flag_periodic_box)
                      r_ij=(REAL)calc_sep_periodic((double)(x_particles_local[cell_index_sort[i_particle]]),
                                                   (double)(y_particles_local[cell_index_sort[i_particle]]),
                                                   (double)(z_particles_local[cell_index_sort[i_particle]]),
                                                   (double)(x_particles_buffer[j_particle]),
                                                   (double)(y_particles_buffer[j_particle]),
                                                   (double)(z_particles_buffer[j_particle]),
                                                   (double)box_size);
                    else
                      r_ij=(REAL)add_quad(3,
                                          (double)(x_particles_local[cell_index_sort[i_particle]]-x_particles_buffer[j_particle]),
                                          (double)(y_particles_local[cell_index_sort[i_particle]]-y_particles_buffer[j_particle]),
                                          (double)(z_particles_local[cell_index_sort[i_particle]]-z_particles_buffer[j_particle]));
                    // Insert into neighbour list (if need be)
                    if(n_neighbours[k_particle]>0){
                      i_neighbour=MIN(N_neighbours,n_neighbours[k_particle]); // Initially give it a place at end of the list
                      while(r_ij<r_neighbours[k_particle][i_neighbour-1]){
                        // Shift particles up if r_ij is smaller than 
                        //   a radius already in the list.  Don't shift
                        //   the N_neighbours-1'th item (it will be overwritten instead)
                        if(i_neighbour<N_neighbours){
                          r_neighbours[k_particle][i_neighbour]  =r_neighbours[k_particle][i_neighbour-1];
                          if(flag_use_velocities)
                            v_neighbours[k_particle][i_neighbour]=v_neighbours[k_particle][i_neighbour-1];
                          if(flag_use_mass && flag_multimass)
                            m_neighbours[k_particle][i_neighbour]=m_neighbours[k_particle][i_neighbour-1];
                        }
                        i_neighbour--; // **
                        if(i_neighbour==0) break;
                      }
                    }
                    else
                      i_neighbour=0;
                    // This test succeeds only if '**' is reached or
                    //   if n_neighbours[k_particle]<N_neighbours
                    if(i_neighbour<N_neighbours){
                      r_neighbours[k_particle][i_neighbour] =r_ij;
                      if(flag_use_velocities)
                        v_neighbours[k_particle][i_neighbour]=(REAL)add_quad(3,
                                                                             (double)(vx_particles_buffer[j_particle]),
                                                                             (double)(vy_particles_buffer[j_particle]),
                                                                             (double)(vz_particles_buffer[j_particle]));
                      if(flag_use_mass && flag_multimass)
                        m_neighbours[k_particle][i_neighbour]=m_particles_buffer[j_particle];
                      n_neighbours[k_particle]++;
                    }
                  } // scan over exchanged particles
                } // if status>0
              } // scan over cell's particles
              // Once all the cells in the buffer have been scanned, check the status of each particle in the cell
              for(i_particle=start_particle,k_particle=0;k_particle<n_p_cell_local;i_particle++,k_particle++){
                if(n_neighbours[k_particle]>=N_neighbours && status[k_particle]>0){
                  // This must be the first buffer scan that checked >=N_neighbour particles.  To be sure we found all the
                  //   neighbours, we must scan at least sqrt(3) times this radius (because scanning a box favours diagonals)
                  if(status[k_particle]>2*grid_size)
                    status[k_particle]=MAX(2,buffer_size);
                  else
                    status[k_particle]--;
                }
                // All neighbours should have been found if status==0 ... perform final calculations
                if(status[k_particle]==0){
                  r_smooth[i_particle]=r_neighbours[k_particle][N_neighbours-1];
                  if(flag_calc_sigma_v)
                    calc_stat(v_neighbours[k_particle],NULL,N_neighbours,ADaPS_REAL,CALC_STAT_STDDEV,&(sigma_v_smooth[i_particle]));
                  if(flag_calc_rho){
                    rho_temp=0.;
                    inv_r_smooth3=pow(r_smooth[i_particle],-3.);
                    if(flag_multimass){
                      for(l_particle=0;l_particle<N_neighbours;l_particle++){
                    	  if(r_neighbours[k_particle][l_particle]<r_smooth[i_particle]){
                          f_table =r_neighbours[k_particle][l_particle]/r_smooth[i_particle];
                          i_table =(int)(f_table*(double)N_KERNEL_TABLE);
                          w_kernel=inv_r_smooth3*(kernel_table[i_table]+(kernel_table[i_table+1]-kernel_table[i_table])*(f_table-kernel_radius[i_table])*N_KERNEL_TABLE);
                    	    rho_temp+=m_neighbours[k_particle][l_particle]*w_kernel; 
                        }      
                    	}
                    }
                    else{
                      for(l_particle=0;l_particle<N_neighbours;l_particle++){
                    	  if(r_neighbours[k_particle][l_particle]<r_smooth[i_particle]){
                          f_table =r_neighbours[k_particle][l_particle]/r_smooth[i_particle];
                          i_table =(int)(f_table*(double)N_KERNEL_TABLE);
                          w_kernel=inv_r_smooth3*(kernel_table[i_table]+(kernel_table[i_table+1]-kernel_table[i_table])*(f_table-kernel_radius[i_table])*N_KERNEL_TABLE);
                    	    rho_temp+=w_kernel; 
                        }      
                    	}
                      rho_temp*=m_p;
                    }
                    rho_smooth[i_particle]=(REAL)rho_temp;
                  }
                  status[k_particle]--;
                  n_left_cell_local--;
                  n_left_total_local--;
                }
              } // check status 
            } // loop over i_rank
            if(n_left_total<n_next_update){
              time(&time_elapsed);time_elapsed-=time_start;
              time_remaining=((double)(n_particles)/(double)(n_particles-n_left_total)-1.)*time_elapsed;
              seconds2ascii(time_remaining,time_string);
              if(SID.n_proc>1)
                SID_log("%3d%% (%lld particles and %d cores) remaining (ETR=%s)",
                        SID_LOG_COMMENT,
                        (int)(100.*(double)n_left_total/(double)n_particles),
                        n_left_total,
                        n_ranks_left,
                        time_string);
              else
                SID_log("%3d%% (%lld particles) remaining (ETR=%s)",
                        SID_LOG_COMMENT,
                        (int)(100.*(double)n_left_total/(double)n_particles),
                        n_left_total,
                        time_string);
              n_next_update-=n_update_step;
            }
          } // while n_left_cell_local>0
        } // loop over cells in z
      } // loop over cells in y
    } // loop over cells in x
    flag_rank_alive=0;
#ifdef USE_MPI
    MPI_Allreduce(&flag_rank_alive,   &n_ranks_left,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&n_left_total_local,&n_left_total,1,MPI_SIZE_T, MPI_SUM,MPI_COMM_WORLD);
#else
    n_ranks_left=flag_rank_alive;
    n_left_total=n_left_total_local;
#endif
    // Continue exchanging so long as there is still an active rank (even if this rank is done)
    if(n_ranks_left>0){
      for(i_rank=0;i_rank<SID.n_proc;i_rank++){
        exchange_buffer(i_rank,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        n_grid,
                        flag_periodic_box,
                        flag_use_velocities,
                        flag_use_mass,
                        flag_multimass,
                        cell_index_local,
                        cell_index_sort,
                        n_particles_local,
                        x_particles_local,
                        y_particles_local,
                        z_particles_local,
                        vx_particles_local,
                        vy_particles_local,
                        vz_particles_local,
                        m_particles_local,
                        x_particles_buffer,
                        y_particles_buffer,
                        z_particles_buffer,
                        vx_particles_buffer,
                        vy_particles_buffer,
                        vz_particles_buffer,
                        m_particles_buffer,
                        &n_particles_buffer,
                        n_buffer_max);
      }
      if(n_left_total<n_next_update){
        time(&time_elapsed);time_elapsed-=time_start;
        time_remaining=((double)(n_particles)/(double)(n_particles-n_left_total)-1.)*time_elapsed;
        seconds2ascii(time_remaining,time_string);
        if(SID.n_proc>1)
          SID_log("%3d%% (%lld particles and %d cores) remaining (ETR=%s)",
                  SID_LOG_COMMENT,
                  (int)(100.*(double)n_left_total/(double)n_particles),
                  n_left_total,
                  n_ranks_left,
                  time_string);
        else
          SID_log("%3d%% (%lld particles) remaining (ETR=%s)",
                  SID_LOG_COMMENT,
                  (int)(100.*(double)n_left_total/(double)n_particles),
                  n_left_total,
                  time_string);
        n_next_update-=n_update_step;
      }
    }
  } // loop until n_ranks_left==0

  // Free the neighbour-info arrays
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  for(j_particle=0;j_particle<n_p_cell_max;j_particle++){
    SID_free((void **)&r_neighbours[j_particle]);
    if(flag_use_velocities)
      SID_free((void **)&v_neighbours[j_particle]);
    if(flag_use_mass && flag_multimass)
      SID_free((void **)&m_neighbours[j_particle]);
  }
  SID_free((void **)&r_neighbours);
  if(flag_use_velocities)
    SID_free((void **)&v_neighbours);
  if(flag_use_mass && flag_multimass)
    SID_free((void **)&m_neighbours);
  SID_free((void **)&n_neighbours);
  SID_free((void **)&status);

  // Free exchange buffer arrays
  SID_free((void **)&x_particles_buffer);
  SID_free((void **)&y_particles_buffer);
  SID_free((void **)&z_particles_buffer);
  if(flag_use_velocities){
    SID_free((void **)&vx_particles_buffer);
    SID_free((void **)&vy_particles_buffer);
    SID_free((void **)&vz_particles_buffer);
  }
  if(flag_use_mass && flag_multimass)
    SID_free((void **)&m_particles_buffer);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Storing arrays...",SID_LOG_OPEN);
  ADaPS_store(&(plist->data),(void *)r_smooth,"r_smooth_%s",ADaPS_DEFAULT,species_name);
  if(flag_calc_sigma_v)
    ADaPS_store(&(plist->data),(void *)sigma_v_smooth,"sigma_v_%s",ADaPS_DEFAULT,species_name);
  if(flag_calc_rho)
    ADaPS_store(&(plist->data),(void *)rho_smooth,"rho_%s",ADaPS_DEFAULT,species_name);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Done.",SID_LOG_CLOSE);
}



            

