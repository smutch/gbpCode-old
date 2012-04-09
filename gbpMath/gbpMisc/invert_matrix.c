#include <string.h>
#include <gbpLib.h>
#include <gbpMisc.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

void invert_square_matrix(double *matrix_in,
                          int     size,
                          double *inverse_out){
   gsl_matrix       *matrix; 
   gsl_matrix       *LU_decomp; 
   gsl_matrix       *inverse; 
   gsl_permutation  *permutation;
   int               signum;
   gsl_matrix_view   m;
   double           *matrix_copy;
   int i_x,i_y,i_M;

   // Allocate some RAM and set some pointers
   matrix_copy=(double *)SID_malloc(sizeof(double)*size*size);
   memcpy(matrix_copy,matrix_in,sizeof(double)*size*size);
   m          =gsl_matrix_view_array(matrix_copy,size,size);
   matrix     =&m.matrix;
   LU_decomp  =matrix;
   inverse    =gsl_matrix_alloc(size,size);
   permutation=gsl_permutation_alloc(size);

   // Perform LU decomposition.  Result ends-up in matrix, 
   //   which is pointed to by LU_decomp
   gsl_linalg_LU_decomp(matrix,   permutation,&signum);
   gsl_linalg_LU_invert(LU_decomp,permutation,inverse);

   // Copy result into the out-going array
   for(i_y=0,i_M=0;i_y<size;i_y++){
      for(i_x=0;i_x<size;i_x++,i_M++){
         inverse_out[i_M]=gsl_matrix_get(inverse,i_x,i_y);
      }
   }

   // Test that the matrix*inverse=identity
   /*
   for(i_x=0;i_x<size;i_x++){
      for(i_y=0;i_y<size;i_y++,i_M++){
         double test;
         int j_M;
         test=0.;
         for(j_M=0;j_M<size;j_M++) test+=matrix_in[i_y*size+j_M]*inverse_out[j_M*size+i_x];
         fprintf(stderr,"%le ",test);
      }
      fprintf(stderr,"\n");
   }
   */

   // Free some RAM
   gsl_matrix_free(inverse);
   gsl_permutation_free(permutation);
   SID_free(SID_FARG matrix_copy);
}

