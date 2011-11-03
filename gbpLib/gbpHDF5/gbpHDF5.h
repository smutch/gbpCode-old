#ifndef GBPHDF5_AWAKE
#define GBPHDF5_AWAKE
#include <hdf5.h>

// Function declarations
#ifdef __cplusplus
extern "C" {
#endif
int read_HDF5_1D(char *arrayName, char *filename, float **data_out, float *themin, float *themax, int *count);
int read_HDF5_3D(char *arrayName, char *filename, float **data_out, float *themin, float *themax, int *count);
#ifdef __cplusplus
}
#endif
#endif
