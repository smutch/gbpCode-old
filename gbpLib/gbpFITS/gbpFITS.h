#ifndef GBPFITS_AWAKE
#define GBPFITS_AWAKE
#include <fitsio.h>

int read_image_FITS(void **image,SID_Datatype *dtype,int *n_D,int **D,char *filename);
int write_image_FITS(void *image,SID_Datatype dtype,int n_D,long *D,char *filename);

#endif
