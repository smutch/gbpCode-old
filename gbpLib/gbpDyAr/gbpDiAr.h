#ifndef DiAr_ON
#define DiAr_ON


// Function declarations
#ifdef __cplusplus
extern "C" {
#endif
void DiAr_init(DiAr **array);
void DiAr_free(void **array);

void *DiAr_get(DiAr *array,size_t element_size,size_t element);
void  DiAr_set(DiAr *array,size_t element_size,size_t element,void *ptr);

double DiAr_get_double(DiAr *array,size_t element);
void   DiAr_set_double(DiAr *array,size_t element,double value);


#ifdef __cplusplus
}
#endif
#endif
