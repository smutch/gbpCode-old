#ifndef GBPPHK_AWAKE
#define GBPPHK_AWAKE
#include <math.h>
#include <hilbert_Moore.h>
#include <SK_sfc_curve.h>
#include <SK_sfc_boundary.h>
#include <SK_hilbert_util.h>
typedef sfc_key_t PHK_t;

PHK_t compute_PHK_from_Cartesian(int PHK_bit_size,int n_D,...);

#endif

