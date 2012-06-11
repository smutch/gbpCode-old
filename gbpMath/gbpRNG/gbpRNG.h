#ifndef GBPRNG_AWAKE
#define GBPRNG_AWAKE
#if USE_SPRNG
  #if USE_MPI == 0
    #define  SIMPLE_SPRNG
  #endif
  #include <sprng.h>
#endif

#define RNG_DEFAULT 2 
#define RNG_GLOBAL  4 // Create the same stream for all ranks

typedef struct RNG_info RNG_info;
struct RNG_info{
  double  GaussBak;
  int     IGauss;
  int     seed;
#if USE_SPRNG
  int    *stream;
#else
  long    stream;
#endif
  int     initialized;
  int     global;
};

// Function definitions
void    init_RNG(int *seed,RNG_info *RNG,int mode);
void    free_RNG(RNG_info *RNG);
void    init_seed_from_clock(int *seed);
GBPREAL random_number(RNG_info *RNG);
GBPREAL random_gaussian(RNG_info *RNG);
GBPREAL random_lognormal(RNG_info *RNG,double mu,double sigma);
float   ran1(long *idum);
void add_gaussian_noise(double *data,int n_data,int *seed,double sigma,double *covariance);

#endif
