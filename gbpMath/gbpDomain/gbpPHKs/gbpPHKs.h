#ifndef GBPPHK_AWAKE
#define GBPPHK_AWAKE
#include <math.h>
#include <hilbert_Moore.h>
#include <SK_sfc_curve.h>
#include <SK_sfc_boundary.h>
#include <SK_hilbert_util.h>
typedef sfc_key_t PHK_t;

#define PHK_DIM_SIZE(N) (1L << N)

#define SID_PHK_T SID_SIZE_T

PHK_t compute_PHK_from_Cartesian(int PHK_bit_size,int n_D,...);
void compute_PHK_volume_keys(int PHK_bit_size,PHK_t key_in,int shell_min,int shell_max,int *n_keys_return,PHK_t **keys_return_in);
void compute_PHK_boundary_keys(int PHK_bit_size,PHK_t domain_key_min,PHK_t domain_key_max,int *n_keys_return,PHK_t **keys_return_in);
void compute_PHK_to_Cartesian(int PHK_bit_size,PHK_t key,int *i_x,int *i_y,int *i_z);

#endif

