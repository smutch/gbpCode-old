#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

void init_seed_from_clock(int *seed){
  time_t seed_temp;
  time(&seed_temp);
  (*seed)=(int)seed_temp;
}
