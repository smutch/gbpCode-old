#ifndef GBPMISC_AWAKE
#define GBPMISC_AWAKE
#include <gbpLib.h>
#include <gbpInterpolate.h>

typedef struct bisect_af_params bisect_af_params;
struct bisect_af_params{
  interp_info *interp;
  double       value;
};

#define  CENTROID3D_MODE_STEP            2
#define  CENTROID3D_MODE_FACTOR          4
#define  CENTROID3D_MODE_INPLACE         8
#define  CENTROID3D_MODE_RETURN_INDICES 16

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif

void invert_square_matrix(double *matrix_in,
                          int     size,
                          double *inverse_out);
double take_aln(double val);
double take_alog10(double val);
double take_ln(double val);
double take_log10(double val);
double add_quad(int n_d, ...);
void apply_rotation(double  x_hat,
                    double  y_hat,
                    double  z_hat,
                    double  theta,
                    double *x_i,
                    double *y_i,
                    double *z_i);
void   compute_Daubechies_scaling_fctns(int D_order,int l_max,double **x_return,double **y_return,int *n_return);
void   force_periodic(GBPREAL *coord,GBPREAL min,GBPREAL box_size);
void   force_periodic_double(double *coord,double min,double box_size);
double d_periodic(double d,double box_size);
int  compute_centroid3D(double  *W,
                        double  *x_in,
                        double  *y_in,
                        double  *z_in,
                        int      n,
                        double   R_min,
                        double   step,
                        int      convergence_N,
                        int      mode,
                        double  *xcen_out,
                        double  *ycen_out,
                        double  *zcen_out);
void   compute_triaxiality(double     *x_in,
                           double     *y_in,
                           double     *z_in,
                           double      x_cen,
                           double      y_cen,
                           double      z_cen,
                           double      box_size,
                           int         n_particles,
                           size_t     *sort_index,
                           double      return_values[3],
                           double      return_vectors[3][3]);
size_t find_index(size_t *y,size_t y_find,size_t  n,size_t *index);
int    find_index_int(int *y,int y_find,int  n,size_t *index);
int    is_a_member(void *candidate,void *list,int n_list,SID_Datatype type);
double bisect_array(interp_info *interp,
                    double       value,
                    double       threshold);
void init_array_linear(double   val_min,
                       double   val_max,
                       int      n_val,
                       double **val,
                       double  *step);
void init_array_log(double   val_min,
                    double   val_max,
                    int      n_val,
                    double **val);
#ifdef __cplusplus
}
#endif
#endif

