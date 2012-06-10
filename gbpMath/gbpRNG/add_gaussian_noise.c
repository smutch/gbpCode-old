#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMisc.h>
#include <gbpRNG.h>
#include <gsl/gsl_linalg.h>

void add_gaussian_noise(double *data,int n_data,int *seed,double sigma,double *covariance){
  int         i_data;
  int         j_data;
  double     *data_temp;
  RNG_info   *RNG;

  // Initialize RNG
  RNG=(RNG_info *)SID_malloc(sizeof(RNG_info));
  init_RNG(seed,RNG,RNG_DEFAULT);

  // Generate a random displacement vector
  gsl_vector *b;
  b=gsl_vector_calloc(n_data);
  for(i_data=0;i_data<n_data;i_data++)
    gsl_vector_set(b,i_data,sigma*random_gaussian(RNG));

  // Use the rotated covariance matrix (if available)
  if(covariance!=NULL){
    // Perform Cholesky Decomposition
    gsl_matrix *m;
    m=gsl_matrix_calloc(n_data,n_data);
    for(i_data=0;i_data<n_data;i_data++)
       for(j_data=0;j_data<n_data;j_data++)
          gsl_matrix_set(m,i_data,j_data,covariance[i_data*n_data+j_data]);
    gsl_linalg_cholesky_decomp(m);
    for(i_data=0;i_data<n_data;i_data++)
       for(j_data=i_data+1;j_data<n_data;j_data++)
          gsl_matrix_set(m,i_data,j_data,0.);
    for(i_data=0;i_data<n_data;i_data++){
       for(j_data=0;j_data<n_data;j_data++)
          data[i_data]+=gsl_matrix_get(m,i_data,j_data)*gsl_vector_get(b,j_data);
    }
    gsl_matrix_free(m);
  }
  else{
    for(j_data=0;j_data<n_data;j_data++)
      data[j_data]+=gsl_vector_get(b,j_data);       
  }

  // Clean-up
  free_RNG(RNG);
  SID_free(SID_FARG RNG);
  gsl_vector_free(b);
}

