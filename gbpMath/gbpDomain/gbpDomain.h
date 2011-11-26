#ifndef GBPDOMAIN_AWAKE
#define GBPDOMAIN_AWAKE
#ifndef GBPFFTW_AWAKE
  #define GBPFFTW_AWAKE
  #ifdef USE_MPI
    #ifdef USE_DOUBLE
      #include <drfftw_mpi.h>
    #else
      #include <srfftw_mpi.h>
    #endif
  #else
    #ifdef USE_DOUBLE
      #include <drfftw.h>
    #else
      #include <srfftw.h>
    #endif
  #endif
#endif

#include <gbpPHKs.h>

typedef struct slab_info slab_info;
struct slab_info{
  double x_min_local;
  double x_max_local;
  double x_max;
  int    n_x_local;
  int    i_x_start_local;
  int    i_x_stop_local;
  int    rank_to_left;
  int    rank_to_right;
};

typedef struct field_info field_info;
struct field_info{
  // Array storing the field
  fftw_real        *field_local;
  fftw_complex     *cfield_local;
  // Field domain
  double           *L;
  double           *dR;
  double           *dk;
  double          **R_field;
  double          **k_field;
  double           *k_Nyquist;
  // Field sizes
  int               n_d;
  int              *n;
  int              *n_R_local;
  int              *n_k_local;
  int              *i_R_start_local;
  int              *i_k_start_local;
  int              *i_R_stop_local;
  int              *i_k_stop_local;
  size_t            n_field;
  size_t            n_field_R_local;
  size_t            n_field_k_local;
  size_t            total_local_size;
  int               pad_size_R;
  int               pad_size_k;
  // flags
  int               flag_padded;
  // FFTW plans
#ifdef USE_MPI
  rfftwnd_mpi_plan  plan;
  rfftwnd_mpi_plan  iplan;
#else
  rfftwnd_plan      plan;
  rfftwnd_plan      iplan;
#endif
  // Slab info
  slab_info         slab;
};

void init_field(int       n_d,
                int      *n,
                double   *L,
                field_info *FFT);
void free_field(field_info *FFT);
void clear_field(field_info *FFT);
void set_exchange_ring_ranks(int *rank_to,
                             int *rank_from,
                             int  i_rank);
void exchange_ring_buffer(void     *send_buffer,
                          size_t    buffer_type_size,
                          size_t    send_count,
                          void     *receive_buffer,
                          size_t   *receive_count,
                          int       i_rank);
void exchange_slab_buffer_left(void      *send_buffer,
                               size_t     send_buffer_size,
                               void      *receive_buffer,
                               size_t    *receive_buffer_size,
                               slab_info *slab);
void exchange_slab_buffer_right(void      *send_buffer,
                                size_t     send_buffer_size,
                                void      *receive_buffer,
                                size_t    *receive_buffer_size,
                                slab_info *slab);

#endif

