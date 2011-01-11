#include <gsl/gsl_linalg.h>

#define MCMC_COVERAGE_LINEAR 0
#define MCMC_COVERAGE_LOG    1
#define MCMC_NAME_SIZE       256

#define MCMC_MODE_DEFAULT    0
#define MCMC_MODE_SERIAL     1
#define MCMC_MODE_PARALLEL   2
#define MCMC_REPORT_PROPS    4
#define MCMC_NO_CVM_UPDATE   8

#define MCMC_MAP_RETURN_GOOD 0
#define MCMC_MAP_RETURN_BAD  1

typedef struct MCMC_DS_info MCMC_DS_info;
struct MCMC_DS_info {
  char          name[MCMC_NAME_SIZE];
  int           n_M;
  double       *M_target;
  double       *dM_target;
  double      **array;
  char        **array_name;
  void         *params;
  int           n_arrays;
  MCMC_DS_info *next;
};

typedef struct MCMC_info MCMC_info;
struct MCMC_info {
  char     filename_output_dir[256];
  char    *problem_name;
  int     (*map_P_to_M)(double *,MCMC_info *,double **);
  void    (*compute_MCMC_ln_likelihood)(MCMC_info *MCMC,double **M,double *ln_likeliood);
  void    *params;
  char   **P_names;
  double  *P_init;
  int     *flag_limit;
  double  *P_limit_min;
  double  *P_limit_max;
  int      n_arrays;
  char   **array_name;
  double **array;
  double  *V;
  gsl_matrix *m;
  int      n_P;
  int      n_DS;
  int      n_M_total;
  double   temperature;
  int      n_avg;
  int      n_avg_covariance;
  int      n_iterations;
  int      n_iterations_burn;
  int      n_thin;
  int      coverage_size;
  int      flag_autocor_on;
  int      flag_integrate_on;
  int      flag_analysis_on;
  int      first_map_call;
  int      mode;
  MCMC_DS_info *DS;
  MCMC_DS_info *last;
};

void init_MCMC(MCMC_info *MCMC,const char *problem_name,void *params,int (*f)(double *,MCMC_info *,double **),int n_P,double *P_init,char **P_names,int *flag_limit,double *P_limit_min,double *P_limit_max,int n_arrays,...);
void add_DS_to_MCMC(MCMC_info *MCMC,const char *name,int n_M,double *DS,double *dDS,void *params,int n_arrays,...);
void add_covariance_to_MCMC(MCMC_info *MCMC,double *V);
void free_MCMC(MCMC_info *MCMC);
void generate_new_MCMC_link(MCMC_info *MCMC,double *P_old,int n_P,double temperature,RNG_info *RNG,gsl_matrix *m,gsl_vector *b,double *P_new);
void compute_MCMC_ln_likelihood_default(MCMC_info *MCMC,double **M,double *ln_likeliood);
void compute_MCMC_chain_stats(double **x,int n_x,int n_avg,double *x_min,double *x_bar_in,double *x_max,double *x_sigma,double **auto_cor,double *slopes,double *dP_sub,double *ln_likelihood_in,double *ln_likelihood_min,double *ln_likelihood_avg,double *ln_likelihood_max);
void compute_MCMC(MCMC_info *MCMC);
void set_MCMC_temperature(MCMC_info *MCMC,double temperature);
void set_MCMC_iterations(MCMC_info *MCMC,int n_avg,int n_avg_covariance,int n_thin,int n_burn,int n_integrate);
void set_MCMC_coverage_size(MCMC_info *MCMC,int coverage_size);
void set_MCMC_likeliood_function(MCMC_info *MCMC,void (*likelihood_function)(MCMC_info *,double **,double *));
void set_MCMC_directory(MCMC_info *MCMC,char *directory);
void set_MCMC_mode(MCMC_info *MCMC,int mode);

