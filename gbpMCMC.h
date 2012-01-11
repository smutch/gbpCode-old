#ifndef MCMC_AWAKE
#define MCMC_AWAKE

#include <gbpRNG.h>
#include <gsl/gsl_linalg.h>

#define MCMC_COVERAGE_LINEAR 0
#define MCMC_COVERAGE_LOG    1
#define MCMC_NAME_SIZE       256

#define MCMC_MODE_SERIAL           1
#define MCMC_MODE_PARALLEL         2
#define MCMC_MODE_AUTOTUNE         8
#define MCMC_MODE_NO_MAP_WRITE    16
#define MCMC_MODE_REPORT_PROPS    32
#define MCMC_MODE_REFLECT_PRIORS  64
#define MCMC_MODE_MINIMIZE_IO    128
#define MCMC_MODE_DEFAULT        MCMC_MODE_AUTOTUNE|MCMC_MODE_SERIAL

#define MCMC_MAP_RETURN_GOOD 0
#define MCMC_MAP_RETURN_BAD  1

#define MCMC_DEFAULT_SUCCESS_TARGET          35.
#define MCMC_DEFAULT_SUCCESS_THRESH           5.
#define MCMC_DEFAULT_COVARIANCE_THRESH       10.
#define MCMC_DEFAULT_N_AUTOTUNE               1
#define MCMC_DEFAULT_N_AUTOTUNE_RANDOMIZE     0
#define MCMC_DEFAULT_N_AUTOTUNE_TEMPERATURE 100
#define MCMC_DEFAULT_N_AUTOTUNE_COVMTX      (10*MCMC->n_P*MCMC->n_P*(1+(int)(100./MCMC->success_target)))
#define MCMC_AUTOTUNE_CONVERGENCE_THRESH     1e-6

typedef struct MCMC_DS_info MCMC_DS_info;
struct MCMC_DS_info {
  char          name[MCMC_NAME_SIZE];
  int           n_D;
  int          *D;
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
  SID_Comm *comm;
  char      filename_output_dir[MAX_FILENAME_LENGTH];
  char      problem_name[MCMC_NAME_SIZE];
  int       (*map_P_to_M)(double *,MCMC_info *,double **);
  void      (*compute_MCMC_ln_likelihood)(MCMC_info *MCMC,double **M,double *P,double *ln_likeliood_DS,int *n_DoF_DS,double *ln_likeliood_all,int *n_DoF_all);
  int       my_chain;
  int       n_chains;
  void    *params;
  char   **P_names;
  double  *P_init;
  double  *P_new;
  double  *P_last;
  double  *P_chain;
  double  *P_limit_min;
  double  *P_limit_max;
  double  *P_min;
  double  *P_max;
  double  *P_avg;
  double  *dP_avg;
  double  *P_best;
  double  *P_peak;
  double  *P_lo_68;
  double  *P_hi_68;
  double  *P_lo_95;
  double  *P_hi_95;
  double   ln_Pr_new;
  double   ln_Pr_last;
  double   ln_Pr_chain;
  int     *n_M;
  double **M_new;
  double **M_last;
  int      n_arrays;
  char   **array_name;
  double **array;
  double  *V;
  gsl_matrix *m;
  gsl_vector *b;
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
  int      flag_no_map_write;
  int      flag_integrate_on;
  int      flag_analysis_on;
  int      flag_init_chain;
  int      first_map_call;
  int      first_link_call;
  size_t   n_map_calls;
  int      mode;
  int      n_success;
  int      n_fail;
  int      n_propositions;
  double   success_target;
  double   success_threshold;
  double   covariance_threshold;
  int      n_autotune;   
  int      n_autotune_randomize;   
  int      n_autotune_temperature;   
  int      n_autotune_covariance;
  int      first_parameter_call;
  int      first_chain_call;
  int      first_likelihood_call;
  double   ln_likelihood_last;
  double   ln_likelihood_new;
  double   ln_likelihood_best;
  double   ln_likelihood_peak;
  double   ln_likelihood_chain;
  double  *ln_likelihood_DS;
  double  *ln_likelihood_DS_best;
  double  *ln_likelihood_DS_peak;
  int     *n_DoF_DS;
  int     *n_DoF_DS_best;
  int     *n_DoF_DS_peak;
  int      n_DoF;
  int      n_DoF_best;
  int      n_DoF_peak;
  int      P_name_length;
  char     P_name_format[8];
  int      seed;
  RNG_info *RNG;
  double       *P;
  MCMC_DS_info *DS;
  MCMC_DS_info *last;
  // This stuff is used when MCMC_MODE_MINIMIZE_IO is on
  char   *flag_success_buffer;
  double *ln_likelihood_new_buffer;
  double *P_new_buffer;
  double *M_new_buffer;
};

#ifdef __cplusplus
extern "C" {
#endif
void init_MCMC(MCMC_info  *MCMC,
               const char *problem_name,
               void       *params,
               int         (*f)(double *,MCMC_info *,double **),
               int         n_P,
               double     *P_init,
               char      **P_names,
               double     *P_limit_min,
               double     *P_limit_max,
               int         n_arrays,...);
void start_loop_MCMC(MCMC_info *MCMC);
void end_loop_MCMC(MCMC_info *MCMC);
void restart_MCMC(MCMC_info *MCMC);
void init_MCMC_arrays(MCMC_info *MCMC);
void free_MCMC_DS(MCMC_info *MCMC);
void free_MCMC_arrays(MCMC_info *MCMC);
void free_MCMC_covariance(MCMC_info *MCMC);
void free_MCMC(MCMC_info *MCMC);
void add_MCMC_DS(MCMC_info *MCMC,const char *name,int n_D,int *D,double *DS,double *dDS,void *params,int n_arrays,...);
void autotune_MCMC(MCMC_info *MCMC);
void autotune_MCMC_randomize(MCMC_info *MCMC);
void autotune_MCMC_temperature(MCMC_info *MCMC);
void autotune_MCMC_covariance(MCMC_info *MCMC);
void compute_MCMC_ln_likelihood_default(MCMC_info *MCMC,double **M,double *P,double *ln_likeliood_DS,int *n_DoF_DS,double *ln_likeliood_all,int *n_DoF_all);
void compute_MCMC_chain_stats(double **x,
                              int      n_x,
                              int      n_avg,
                              double  *x_min,
                              double  *x_bar_in,
                              double  *x_max,
                              double  *x_sigma,
                              double **auto_cor,
                              double  *slopes,
                              double  *dP_sub,
                              double  *ln_likelihood_in,
                              double  *ln_likelihood_min,
                              double  *ln_likelihood_avg,
                              double  *ln_likelihood_max);
void compute_MCMC(MCMC_info *MCMC);
void analyze_MCMC(MCMC_info *MCMC);
void set_MCMC_covariance(MCMC_info *MCMC,double *V);
void set_MCMC_temperature(MCMC_info *MCMC,double temperature);
void set_MCMC_iterations(MCMC_info *MCMC,int n_avg,int n_thin,int n_burn,int n_integrate);
void set_MCMC_coverage_size(MCMC_info *MCMC,int coverage_size);
void set_MCMC_likelihood_function(MCMC_info *MCMC,void (*likelihood_function)(MCMC_info *,double **,double *,double *,int *,double *,int *));
void set_MCMC_directory(MCMC_info *MCMC,char *directory);
void set_MCMC_mode(MCMC_info *MCMC,int mode);
void set_MCMC_autotune(MCMC_info *MCMC,
                       double     success_target,
                       double     success_threshold,
                       double     covariance_threshold,
                       int        n_autotune,
                       int        n_autotune_randomize,
                       int        n_autotune_temperature,
                       int        n_autotune_covariance);
void read_MCMC_covariance(MCMC_info *MCMC,char *filename);
void read_MCMC_state(MCMC_info *MCMC);
void rm_MCMC_directory(MCMC_info *MCMC);
int  generate_MCMC_chain(MCMC_info *MCMC);
void generate_MCMC_proposition(MCMC_info *MCMC,int flag_chain_init);
void generate_MCMC_parameters(MCMC_info *MCMC);
#ifdef __cplusplus
}
#endif

#endif

