// Free memory allocated by initialize_field()
#include <gbpLib.h>
#include <gbpSort.h>
#include <gbpPHKs.h>

void compute_PHK_boundary_keys(int PHK_bit_size,PHK_t domain_key_min,PHK_t domain_key_max,int *n_keys_return,PHK_t **keys_return){
   hikey_t  basecoord[3];
   hikey_t  coord[3];
   hikey_t  num1dim;
   int      i, j, k;
   int      count;
   unsigned bits;

   // Some type conversions
   bits=(unsigned)PHK_bit_size;

   // Calculate the number of possible Key values per dimension 
   num1dim = PHK_DIM_SIZE(bits);

   // Sanity checks
   if(domain_key_min<0)
     SID_trap_error("domain_key_min<0 (ie %zd<0) in compute_PHK_boundry_keys",ERROR_LOGIC,domain_key_min);
   if(domain_key_max>(num1dim*num1dim*num1dim))
     SID_trap_error("domain_key_max>num1dim^3 (ie %zd>%zd) in compute_PHK_boundry_keys",ERROR_LOGIC,domain_key_max,num1dim*num1dim*num1dim);
   if(domain_key_max<domain_key_min)
     SID_trap_error("domain_key_max<domain_key_min (ie %zd<%zd) in compute_PHK_volume_keys",ERROR_LOGIC,domain_key_min,domain_key_max);

   // Compute the number of keys to be returned
   PHK_t  key;
   PHK_t *PHK_test=NULL;
   int    i_test;
   int    n_test;
   int    flag;
   (*n_keys_return)=0;
   for(key=domain_key_min;key<=domain_key_max;key++){
      compute_PHK_volume_keys(PHK_bit_size,key,0,1,&n_test,&PHK_test);
      for(i_test=0,flag=TRUE;i_test<n_test && flag;i_test++){
         if(PHK_test[i_test]<domain_key_min || PHK_test[i_test]>domain_key_max){
            (*n_keys_return)++;
            flag=FALSE;
         }
      }
   }

   // Always allocate for the key array
   (*keys_return)=(PHK_t *)SID_malloc((*n_keys_return)*sizeof(PHK_t));

   // Find keys
   for(key=domain_key_min,count=0;key<=domain_key_max;key++){
      compute_PHK_volume_keys(PHK_bit_size,key,0,1,&n_test,&PHK_test);
      for(i_test=0,flag=TRUE;i_test<n_test && flag;i_test++){
         if(PHK_test[i_test]<domain_key_min || PHK_test[i_test]>domain_key_max){
            (*keys_return)[count++]=key;
            flag=FALSE;
         }
      }
   }
   SID_free(SID_FARG PHK_test);

   // Sort keys in place
   merge_sort((*keys_return),(*n_keys_return),NULL,SID_PHK_T,SORT_INPLACE_ONLY,SORT_COMPUTE_INPLACE);

}

