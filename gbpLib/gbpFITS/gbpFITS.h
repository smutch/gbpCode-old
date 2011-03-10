#ifndef GBPFITS_AWAKE
#define GBPFITS_AWAKE
#include <fitsio.h>

int read_image_FITS(void **image,SID_Datatype *dtype,int *n_D,int **D,char *filename);
int write_image_FITS(void *image,SID_Datatype dtype,int n_D,int *D,char *filename);
void transpose_array(void *data,SID_Datatype dtype,int n_d,int *d_i);

#endif
