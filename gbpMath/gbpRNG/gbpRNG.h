#ifndef USE_MPI
  #define  SIMPLE_SPRNG
#endif
#include <sprng.h>

#define RNG_DEFAULT 2 
#define RNG_GLOBAL  4 // Create the same stream for all ranks

typedef struct RNG_info RNG_info;
struct RNG_info{
  double  GaussBak;
  int     IGauss;
  int     seed;
  int    *stream;
  int     initialized;
  int     global;
};

// Function definitions
void init_RNG(int *seed,RNG_info *RNG,int mode);
void free_RNG(RNG_info *RNG);
void init_seed_from_clock(int *seed);
REAL random_number(RNG_info *RNG);
REAL random_gaussian(RNG_info *RNG);
REAL random_lognormal(RNG_info *RNG,double mu,double sigma);

