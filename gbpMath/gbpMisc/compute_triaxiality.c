#include <gbpLib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_spline.h>

void compute_triaxiality(double     *x_in,
                         double     *y_in,
                         double     *z_in,
                         double      x_cen,
                         double      y_cen,
                         double      z_cen,
                         double      box_size,
                         int         n_particles,
                         size_t     *sort_index,
                         double      return_values[3],
                         double      return_vectors[3][3]){
  double   s,s_new;
  double   q,q_new;
  double   a_new;
  double   b_new;
  double   c_new;
  int      i,j;
  double   M_tmp;
  double   x_tmp;
  double   y_tmp;
  double   z_tmp;
  double  *m_p;
  double  *x;
  double  *y;
  double  *z;
  double   inv_q2;
  double   inv_s2;
  double   inv_r2,r2;
  double   M[9];
  int      n_iterations;
  int      continue_flag;
  double   convergence     =0.0001;
  int      n_iterations_max=200;
  gsl_vector                *eigen_vector; 
  gsl_matrix                *eigen_matrix; 
  gsl_eigen_symmv_workspace *w; 
  gsl_matrix_view            m; 

  // Initialize a bunch of stuff
  a_new=1.;
  b_new=1.;
  c_new=1.;
  q_new=1.;
  s_new=1.;
  return_vectors[0][0]=1.;return_vectors[1][0]=0.;return_vectors[2][0]=0.;
  return_vectors[0][1]=0.;return_vectors[1][1]=1.;return_vectors[2][1]=0.;
  return_vectors[0][2]=0.;return_vectors[1][2]=0.;return_vectors[2][2]=1.;

  if(n_particles>0)
    continue_flag=TRUE;
  else
    continue_flag=FALSE;
  eigen_vector=gsl_vector_alloc(3); 
  eigen_matrix=gsl_matrix_alloc(3,3); 
  w           =gsl_eigen_symmv_alloc(3); 

  m_p=(double *)SID_malloc(sizeof(double)*n_particles);
  x  =(double *)SID_malloc(sizeof(double)*n_particles);
  y  =(double *)SID_malloc(sizeof(double)*n_particles);
  z  =(double *)SID_malloc(sizeof(double)*n_particles);
  if(sort_index!=NULL){
    for(i=0;i<n_particles;i++){
      m_p[i]=1.;
      x[i]  =d_periodic((double)x_in[sort_index[i]]-x_cen,box_size);
      y[i]  =d_periodic((double)y_in[sort_index[i]]-y_cen,box_size);
      z[i]  =d_periodic((double)z_in[sort_index[i]]-z_cen,box_size);
    }
  }
  else{
    for(i=0;i<n_particles;i++){
      m_p[i]=1.;
      x[i]  =d_periodic((double)x_in[i]-x_cen,box_size);
      y[i]  =d_periodic((double)y_in[i]-y_cen,box_size);
      z[i]  =d_periodic((double)z_in[i]-z_cen,box_size);
    }    
  }

  // Iterate until convergence
  n_iterations =0;
  while(continue_flag){
    q     =q_new;
    s     =s_new;
    inv_q2=1./(q*q);
    inv_s2=1./(s*s);

    // Construct the moment of inertia tensor
    for(i=0;i<9;i++)
      M[i]=0.;
    for(i=0;i<n_particles;i++){
      x_tmp=x[i];
      y_tmp=y[i];
      z_tmp=z[i];
      M_tmp=m_p[i];
      r2   =x_tmp*x_tmp+y_tmp*y_tmp*inv_q2+z_tmp*z_tmp*inv_s2;
      if(r2>0.){
        inv_r2=1./r2;
        M[0]+=x_tmp*x_tmp*M_tmp*inv_r2;
        M[1]+=y_tmp*x_tmp*M_tmp*inv_r2;
        M[2]+=z_tmp*x_tmp*M_tmp*inv_r2;
        M[3]+=x_tmp*y_tmp*M_tmp*inv_r2;
        M[4]+=y_tmp*y_tmp*M_tmp*inv_r2;
        M[5]+=z_tmp*y_tmp*M_tmp*inv_r2;
        M[6]+=x_tmp*z_tmp*M_tmp*inv_r2;
        M[7]+=y_tmp*z_tmp*M_tmp*inv_r2;
        M[8]+=z_tmp*z_tmp*M_tmp*inv_r2;
      }
    }

    // Solve for (and sort) the eigen values and eigen vectors
    m=gsl_matrix_view_array(M,3,3);
    gsl_eigen_symmv(&m.matrix,eigen_vector,eigen_matrix,w); 
    gsl_eigen_symmv_sort(eigen_vector,eigen_matrix,GSL_EIGEN_SORT_ABS_DESC);

    // Convert goofy gsl vectors and such into something simpler to use
    for(i=0;i<3;i++){
      return_values[i]    =sqrt(gsl_vector_get(eigen_vector,i));
      return_vectors[i][0]=gsl_matrix_get(eigen_matrix,i,0);
      return_vectors[i][1]=gsl_matrix_get(eigen_matrix,i,1);
      return_vectors[i][2]=gsl_matrix_get(eigen_matrix,i,2);
    }
    q_new=return_values[1]/return_values[0];
    s_new=return_values[2]/return_values[0];

    // Check for convergence
    n_iterations++;
    if(n_iterations>=n_iterations_max)
      continue_flag=FALSE;
    if((double)fabs((float)((q_new-q)/q))<convergence && (double)fabs((float)((s_new-s)/s))<convergence)
      continue_flag=FALSE;
  }

  // Clean-up
  gsl_eigen_symmv_free(w);
  gsl_vector_free(eigen_vector); 
  gsl_matrix_free(eigen_matrix); 
  SID_free((void **)(&x));
  SID_free((void **)(&y));
  SID_free((void **)(&z));
  SID_free((void **)(&m_p));
}

