
// ADaPS: Adaptive Data Passing Structure
#ifndef ADaPS_ON
#define ADaPS_ON

#include <gbpCommon.h>
#include <gbpSID.h>

#define ADaPS_NAME_LENGTH    40

#define ADaPS_DEFAULT                 2
#define ADaPS_COPY                    4
#define ADaPS_SCALAR_DOUBLE           8
#define ADaPS_SCALAR_FLOAT           16
#define ADaPS_SCALAR_REAL            32
#define ADaPS_SCALAR_SIZE_T          64
#define ADaPS_SCALAR_INT            128
#define ADaPS_COPY_SUBARRAY_DOUBLE  256
#define ADaPS_COPY_SUBARRAY_FLOAT   512
#define ADaPS_COPY_SUBARRAY_REAL   1024
#define ADaPS_COPY_SUBARRAY_SIZE_T 2048
#define ADaPS_COPY_SUBARRAY_INT    4096

#define ADaPS_DOUBLE       0
#define ADaPS_LONG         1
#define ADaPS_FLOAT        2
#define ADaPS_INT          3
#define ADaPS_UINT         4
#define ADaPS_LONG_LONG    5
#define ADaPS_SIZE_T       6
#define ADaPS_STRING       7
#define ADaPS_CUSTOM       8
#if USE_DOUBLE
  #define ADaPS_REAL ADaPS_DOUBLE
#else
  #define ADaPS_REAL ADaPS_FLOAT
#endif

// Define the structure
typedef struct ADaPS_struct ADaPS;
struct ADaPS_struct{
  SID_Datatype data_type;
  char         name[ADaPS_NAME_LENGTH];
  void        *data;
  int          mode;
  ADaPS       *next;
  size_t       data_size;
  void        (*free_function)(void **,void *);
  void        *free_function_params;
};

// Function declarations
#ifdef __cplusplus
extern "C" {
#endif
void ADaPS_init(ADaPS **list);
void ADaPS_store(ADaPS **list,
                 void   *data,
                 const char   *name,
                 int     mode, ...);
void *ADaPS_fetch(ADaPS *list,
                  const char  *name_in,...);
int  ADaPS_exist(ADaPS *list,
                 const char  *name, ...);
void ADaPS_free(void **list);
void ADaPS_deallocate(ADaPS **remove);
void ADaPS_remove(ADaPS **list, 
                  const char   *name,...);
void ADaPS_status(ADaPS *list);
#ifdef __cplusplus
}
#endif
#endif
