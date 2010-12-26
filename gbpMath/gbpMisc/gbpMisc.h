double take_aln(double val);
double take_alog10(double val);
double take_ln(double val);
double take_log10(double val);
void   compute_Daubechies_scaling_fctns(int D_order,int l_max,double **x_return,double **y_return,int *n_return);
void   force_periodic(REAL *coord,REAL min,REAL box_size);
double d_periodic(double d,double box_size);
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

