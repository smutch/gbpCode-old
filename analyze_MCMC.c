#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <gbpLib.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void analyze_MCMC(MCMC_info *MCMC){
  char      filename_output_dir[MAX_FILENAME_LENGTH];
  char      filename_chain_dir[MAX_FILENAME_LENGTH];
  char      filename_results_dir[MAX_FILENAME_LENGTH];
  char      filename_plots_dir[MAX_FILENAME_LENGTH];
  char      filename_run[MAX_FILENAME_LENGTH];
  char      filename_chain[MAX_FILENAME_LENGTH];
  char      filename_chain_iterations[MAX_FILENAME_LENGTH];
  char      filename_stats[MAX_FILENAME_LENGTH];
  char      filename_coverage[MAX_FILENAME_LENGTH];
  char      filename_chain_covariance[MAX_FILENAME_LENGTH];
  char      filename_covariance[MAX_FILENAME_LENGTH];
  char      filename_histograms[MAX_FILENAME_LENGTH];
  char      filename_results[MAX_FILENAME_LENGTH];
  char      filename_stop[MAX_FILENAME_LENGTH];
  char      column_txt[MAX_FILENAME_LENGTH];
  char      problem_name_test[MCMC_NAME_SIZE];
  char      P_name_test[MCMC_NAME_SIZE];
  char      name_test[MCMC_NAME_SIZE];
  char      array_name_test[MCMC_NAME_SIZE];
  char      array_name_temp[MCMC_NAME_SIZE];
  int       n_avg_test;
  int       flag_autocor_on_test;
  int       flag_no_map_write;
  int       n_P_test;
  int       n_DS_test;
  int       n_arrays_test;
  int       n_M_test;
  double    P_init_test;
  double    M_target_test;
  double    dM_target_test;
  int       i_P,j_P,k_P;
  int       i_M,j_M;
  int       ii_P,jj_P;
  int       i_avg;
  int       i_C;
  int       i_DS,n_DS;
  int       i_iteration,j_iteration;
  int       flag_continue;
  char      flag_success;
  int       flag_stop=FALSE;
  int       i_coverage;
  int       j_coverage;
  int       i_proposal;
  int       bin_x;
  int       bin_y;
  int       n_avg;
  int       n_avg_burn;
  int       n_avg_integrate;
  int       n_P,n_C;
  int       n_gals;
  size_t    n_success,n_fail;
  size_t    n_success_all,n_fail_all;
  double   *x_P;
  double  **x_M;
  double   *M_target;
  double   *dM_target;
  double   *V;
  double  **C;
  double   *P_min;
  double   *P_max;
  double   *P_new;
  double   *P_last;
  double   *P_avg;
  double   *dP_avg;
  double  **M_best;
  double  **M_best_parameters;
  double  **M_peak_parameters;
  size_t    n_residual;
  double   *residual;
  double  **M_min;
  double  **M_max;
  double  **M_new;
  double  **M_lo_68;
  double  **M_hi_68;
  double  **M_lo_95;
  double  **M_hi_95;
  size_t ***M_histogram;
  double   **M_last;
  double   **M_avg;
  double   **dM_avg;
  size_t  **P_histogram;
  size_t  **coverage_true;
  size_t  **coverage_false;
  size_t  **coverage_keep;
  double  **mean_L_histogram;
  double  **max_L_histogram;
  uint64_t *temp_out;
  double   **max_L_surface;
  double   **mean_L_surface;
  double    L_x,L_y,L_z;
  int       n_x,n_y,n_z;
  int       n_C_x;
  double    ln_likelihood_new;
  double    ln_likelihood_last;
  double    ln_likelihood_best;
  double    ln_likelihood_peak;
  FILE     *fp_run;
  FILE     *fp_chain;
  FILE     *fp_chain_iterations;
  FILE     *fp_stats;
  FILE     *fp_coverage;
  FILE     *fp_chain_covariance;
  FILE     *fp_covariance;
  FILE     *fp_histograms;
  FILE     *fp_results;
  FILE     *fp_stop;
  int       seed=182743;
  int       i_report;
  int       i_iteration_start;
  int       i_iteration_next_report;
  int       n_accepted;
  int       n_active;
  int       n_iterations;
  int       n_iterations_burn;
  int       n_iterations_integrate;
  int       n_iterations_phase;
  int       n_thin;
  int       n_thin_burn;
  int       n_thin_integrate;
  int       i_thin;
  double    temperature;
  int       flag_autocor_on;
  int       n_coverage;
  int       coverage_size;
  int       i_phase;
  int       i_array;
  int      *n_M;
  int       flag_initialized;
  int       n_used;
  gsl_matrix *m;
  gsl_vector *b;
  double      RN;
  MCMC_DS_info  *current_DS;
  MCMC_DS_info  *next_DS;
  double       ***M_arrays;
  double       ***P_arrays;
  int            *n_P_arrays;
  int            *n_M_arrays;
  double          constant;
  size_t         *histogram_index;
  size_t         *coverage_keep_index;
  size_t          P_best_index;
  size_t          M_best_index;
  double          accumulator;
  size_t          P_lo_68_index;
  size_t          P_hi_68_index;
  size_t          P_lo_95_index;
  size_t          P_hi_95_index;
  size_t          M_lo_68_index;
  size_t          M_hi_68_index;
  size_t          M_lo_95_index;
  size_t          M_hi_95_index;
  double         *P_lo_68;
  double         *P_hi_68;
  double         *P_lo_95;
  double         *P_hi_95;
  double         *P_contour_68;
  double         *P_contour_95;
  double         *P_best;
  double         *P_peak;
  int             my_chain;
  int             flag_init;
  int             i_column;
  int             dummy_t;
  double         *P_i_bar_accum;
  double         *P_ij_bar_accum;
  double         *P_i_bar;
  double         *P_ij_bar;
  double         *V_compute;
  int             n_covariance;
  int             i_covariance;
  int             j_covariance;
  int             n_iterations_file_total;
  int             n_iterations_file_burn;
  int             flag_restart=FALSE;
  int   size_temp1;
  int   size_temp2[2];

  SID_log("Performing MCMC analysis...",SID_LOG_OPEN|SID_LOG_TIMER);

  SID_log("Initializing...",SID_LOG_OPEN);

  // Initialize dataset arrays
  init_MCMC_arrays(MCMC);

  // Initialize some constants etc.
  n_P                   =MCMC->n_P;
  n_DS                  =MCMC->n_DS;
  n_avg                 =MCMC->n_avg;
  my_chain              =MCMC->my_chain;
  n_iterations          =MCMC->n_iterations;
  n_iterations_burn     =MCMC->n_iterations_burn;
  n_iterations_integrate=n_iterations-n_iterations_burn;
  n_thin                =MCMC->n_thin;
  coverage_size         =MCMC->coverage_size;
  flag_autocor_on       =MCMC->flag_autocor_on;
  flag_no_map_write     =MCMC->flag_no_map_write;

  // Initialize arrays used for computing the covariance matrix
  P_i_bar_accum =(double *)SID_malloc(sizeof(double)*n_P);
  P_ij_bar_accum=(double *)SID_malloc(sizeof(double)*n_P*n_P);
  P_i_bar       =(double *)SID_malloc(sizeof(double)*n_P);
  P_ij_bar      =(double *)SID_malloc(sizeof(double)*n_P*n_P);
  V_compute     =(double *)SID_malloc(sizeof(double)*n_P*n_P);
  for(i_P=0,k_P=0;i_P<n_P;i_P++){
    P_i_bar_accum[i_P]=0.;
    P_i_bar[i_P]      =0.;
    for(j_P=0;j_P<n_P;j_P++,k_P++){
      P_ij_bar_accum[k_P]=0.;
      P_ij_bar[k_P]      =0.;
      if(i_P==j_P)
        V_compute[k_P]=1.;
      else
        V_compute[k_P]=0.;
    }
  }

  // Initialize chain arrays etc
  n_M       =MCMC->n_M;
  M_new     =MCMC->M_new;
  M_last    =MCMC->M_last;
  n_M_arrays=(int      *)SID_malloc(sizeof(int      )*n_DS);
  M_arrays  =(double ***)SID_malloc(sizeof(double **)*n_DS);
  current_DS=MCMC->DS;
  i_DS      =0;
  while(current_DS!=NULL){
    next_DS         =current_DS->next;
    n_M[i_DS]       =current_DS->n_M;
    n_M_arrays[i_DS]=current_DS->n_arrays;
    M_arrays[i_DS]  =current_DS->array;
    current_DS      =next_DS;
    i_DS++;
  }
  P_new          =MCMC->P_new;
  P_last         =MCMC->P_last;

  // Set directories
  sprintf(filename_output_dir, "%s/",        MCMC->filename_output_dir);
  sprintf(filename_chain_dir,  "%s/chains/", filename_output_dir);
  sprintf(filename_results_dir,"%s/results/",filename_output_dir);
  sprintf(filename_plots_dir,  "%s/plots/",  filename_output_dir);

  // Set filenames
  sprintf(filename_run,             "%s/run.dat",                  filename_output_dir);
  sprintf(filename_chain,           "%s/chain_trace_%06d.dat",     filename_chain_dir,my_chain);
  sprintf(filename_chain_iterations,"%s/chain_iterations_%06d.dat",filename_chain_dir,my_chain);
  sprintf(filename_chain_covariance,"%s/chain_covariance_%06d.dat",filename_chain_dir,my_chain);
  sprintf(filename_stats,           "%s/chain_stats_%06d.dat",     filename_chain_dir,my_chain);
  sprintf(filename_coverage,        "%s/coverage.dat",             filename_results_dir);
  sprintf(filename_histograms,      "%s/histograms.dat",           filename_results_dir);
  sprintf(filename_covariance,      "%s/covariance.dat",           filename_results_dir);
  sprintf(filename_stop,            "%s/stop",                     filename_output_dir);
  
  // Initialize coverage arrays
  for(i_P=0,n_coverage=0;i_P<n_P;i_P++){
      for(j_P=i_P+1;j_P<n_P;j_P++)
          n_coverage++;
  }
  P_min         =MCMC->P_min;
  P_max         =MCMC->P_max;
  P_avg         =MCMC->P_avg;
  dP_avg        =MCMC->dP_avg;
  P_best        =MCMC->P_best;
  P_peak        =MCMC->P_peak;
  P_lo_68       =MCMC->P_lo_68;
  P_hi_68       =MCMC->P_hi_68;
  P_lo_95       =MCMC->P_lo_95;
  P_hi_95       =MCMC->P_hi_95;
  P_contour_68  =(double  *)SID_malloc(sizeof(double)*n_coverage);
  P_contour_95  =(double  *)SID_malloc(sizeof(double)*n_coverage);
  P_histogram   =(size_t **)SID_malloc(sizeof(size_t *)*n_P);
  mean_L_histogram =(double  **)SID_malloc(sizeof(double *)*n_P);
  max_L_histogram  =(double  **)SID_malloc(sizeof(double *)*n_P);
  for(i_P=0;i_P<n_P;i_P++){
    P_histogram[i_P]     =(size_t *)SID_malloc(sizeof(size_t)*coverage_size);
    mean_L_histogram[i_P]=(double *)SID_malloc(sizeof(double)*coverage_size);
    max_L_histogram[i_P] =(double *)SID_malloc(sizeof(double)*coverage_size);
    P_min[i_P] = DBL_MAX;
    P_max[i_P] =-DBL_MAX;
    P_avg[i_P] = 0.;
    dP_avg[i_P]= 0.;
    for(j_P=0;j_P<coverage_size;j_P++){
      P_histogram[i_P][j_P]     =0;
      mean_L_histogram[i_P][j_P]=0.;
      max_L_histogram[i_P][j_P] =-1e10;
    }
  }
  M_best           =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_best_parameters=(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_peak_parameters=(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_min            =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_max            =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_lo_68          =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_hi_68          =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_lo_95          =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_hi_95          =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_avg            =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  dM_avg           =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_histogram      =(size_t ***)SID_malloc(sizeof(size_t **)*n_DS);
  for(i_DS=0,n_residual=0;i_DS<n_DS;i_DS++){
    M_best[i_DS]           =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_best_parameters[i_DS]=(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_peak_parameters[i_DS]=(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_min[i_DS]            =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_max[i_DS]            =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_lo_68[i_DS]          =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_hi_68[i_DS]          =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_lo_95[i_DS]          =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_hi_95[i_DS]          =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_avg[i_DS]            =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    dM_avg[i_DS]           =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_histogram[i_DS]      =(size_t **)SID_malloc(sizeof(size_t *)*n_M[i_DS]);
    for(i_M=0;i_M<n_M[i_DS];i_M++){
      M_histogram[i_DS][i_M]=(size_t *)SID_malloc(sizeof(size_t)*coverage_size);
      M_avg[i_DS][i_M]  =0.;
      dM_avg[i_DS][i_M] =0.;
      M_min[i_DS][i_M]  = DBL_MAX;
      M_max[i_DS][i_M]  =-DBL_MAX;
      for(j_M=0;j_M<coverage_size;j_M++)
        M_histogram[i_DS][i_M][j_M]=0;
    }
    n_residual=MAX(n_residual,n_M[i_DS]);
  }
  residual         =(double   *)SID_malloc(sizeof(double)*n_residual);
  coverage_true    =(size_t  **)SID_malloc(sizeof(size_t  *)*n_coverage);
  coverage_false   =(size_t  **)SID_malloc(sizeof(size_t  *)*n_coverage);
  coverage_keep    =(size_t  **)SID_malloc(sizeof(size_t  *)*n_coverage);
  max_L_surface    =(double  **)SID_malloc(sizeof(double  *)*n_coverage);
  mean_L_surface   =(double  **)SID_malloc(sizeof(double  *)*n_coverage);
  for(i_coverage=0;i_coverage<n_coverage;i_coverage++){
    coverage_true[i_coverage]    =(size_t *)SID_malloc(sizeof(size_t)*coverage_size*coverage_size);
    coverage_false[i_coverage]   =(size_t *)SID_malloc(sizeof(size_t)*coverage_size*coverage_size);
    coverage_keep[i_coverage]    =(size_t *)SID_malloc(sizeof(size_t)*coverage_size*coverage_size);
    max_L_surface[i_coverage]    =(double *)SID_malloc(sizeof(double)*coverage_size*coverage_size);
    mean_L_surface[i_coverage]   =(double *)SID_malloc(sizeof(double)*coverage_size*coverage_size);
    for(j_coverage=0;j_coverage<coverage_size*coverage_size;j_coverage++){
      coverage_true[i_coverage][j_coverage]    =0;
      coverage_false[i_coverage][j_coverage]   =0;
      coverage_keep[i_coverage][j_coverage]    =0;
      max_L_surface[i_coverage][j_coverage]    =-1e10;
      mean_L_surface[i_coverage][j_coverage]   =0.;
    }
  }
  temp_out=(uint64_t *)SID_malloc(sizeof(uint64_t)*coverage_size*coverage_size);
  SID_log("Done.",SID_LOG_CLOSE);

  // Process the chain(s)
  SID_log("Process chain file(s)...",SID_LOG_OPEN);
  if(flag_no_map_write)
    SID_log("Mapping results are *not* available.",SID_LOG_COMMENT);
  else
    SID_log("Mapping results are available.",SID_LOG_COMMENT);
  if(my_chain==SID.My_rank){
      if((fp_chain=fopen(filename_chain,"rb"))==NULL)
        SID_trap_error("Could not open chain file {%s}.",ERROR_IO_OPEN,filename_chain);
      for(i_iteration=0,flag_init=TRUE,n_used=0;i_iteration<n_iterations;i_iteration++){
        for(i_avg=0;i_avg<n_avg;i_avg++){
          fread(&flag_success,     sizeof(char),  1,  fp_chain);
          fread(&ln_likelihood_new,sizeof(double),1,  fp_chain);
          fread(P_new,             sizeof(double),n_P,fp_chain);
          if(!flag_no_map_write){
            for(i_DS=0;i_DS<n_DS;i_DS++)
              fread(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
          }
          if(i_iteration>=n_iterations_burn){
            if(flag_init || flag_success){
              memcpy(P_last,P_new,(size_t)n_P*sizeof(double));
              ln_likelihood_last=ln_likelihood_new;
              if(!flag_no_map_write){
                for(i_DS=0;i_DS<n_DS;i_DS++)
                  memcpy(M_last[i_DS],M_new[i_DS],(size_t)n_M[i_DS]*sizeof(double));
              }
              flag_init=FALSE;
            }
            n_used++;
            // Compute covariance matrix
            for(i_P=0,k_P=0;i_P<n_P;i_P++){
              P_i_bar_accum[i_P]+=P_last[i_P];
              for(j_P=0;j_P<n_P;j_P++,k_P++)
                P_ij_bar_accum[k_P]+=P_last[i_P]*P_last[j_P];
            }

            // Compute parameter extrema and averages
            for(i_P=0;i_P<n_P;i_P++){
              if(P_last[i_P]<P_min[i_P]) P_min[i_P]=P_last[i_P];
              if(P_last[i_P]>P_max[i_P]) P_max[i_P]=P_last[i_P];
              P_avg[i_P]+=P_last[i_P];
            }

            // Compute mapping-space extrema and averages
            if(!flag_no_map_write){
              for(i_DS=0;i_DS<n_DS;i_DS++){
                for(i_M=0;i_M<n_M[i_DS];i_M++){
                  if(M_last[i_DS][i_M]<M_min[i_DS][i_M]) M_min[i_DS][i_M]=M_last[i_DS][i_M];
                  if(M_last[i_DS][i_M]>M_max[i_DS][i_M]) M_max[i_DS][i_M]=M_last[i_DS][i_M];
                  M_avg[i_DS][i_M]+=M_last[i_DS][i_M];
                }
              }
            }
          }
        }
      }
      // Finish averages
      for(i_P=0;i_P<n_P;i_P++){
        P_avg[i_P]/=(double)n_used;
      }
      if(!flag_no_map_write){
        for(i_DS=0;i_DS<n_DS;i_DS++){
          for(i_M=0;i_M<n_M[i_DS];i_M++)
            M_avg[i_DS][i_M]/=(double)n_used;
        }
      }
      rewind(fp_chain);
      for(i_iteration=0,flag_init=TRUE;i_iteration<n_iterations;i_iteration++){
        for(i_avg=0;i_avg<n_avg;i_avg++){
          fread(&flag_success,     sizeof(char),  1,  fp_chain);
          fread(&ln_likelihood_new,sizeof(double),1,  fp_chain);
          fread(P_new,             sizeof(double),n_P,fp_chain);
          if(!flag_no_map_write){
            for(i_DS=0;i_DS<n_DS;i_DS++)
              fread(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
          }
          if(i_iteration>=n_iterations_burn){
            if(flag_init || flag_success){
              memcpy(P_last,P_new,(size_t)n_P*sizeof(double));      
              ln_likelihood_last=ln_likelihood_new;
              for(i_DS=0;i_DS<n_DS;i_DS++)
                memcpy(M_last[i_DS],M_new[i_DS],(size_t)n_M[i_DS]*sizeof(double));
              flag_init=FALSE;
            }
            // Build coverage maps
            if(flag_success){
               for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
                  bin_x=(double)(coverage_size)*(P_last[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
                  for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
                     bin_y=(double)(coverage_size)*(P_last[j_P]-P_min[j_P])/(P_max[j_P]-P_min[j_P]);
                     if(bin_x>=0 && bin_x<coverage_size && bin_y>=0 && bin_y<coverage_size)
                        coverage_true[i_coverage][bin_x*coverage_size+bin_y]++;
                  }
               }
            }
            else{
               for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
                  bin_x=(double)(coverage_size)*(P_new[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
                  for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
                     bin_y=(double)(coverage_size)*(P_new[j_P]-P_min[j_P])/(P_max[j_P]-P_min[j_P]);
                     if(bin_x>=0 && bin_x<coverage_size && bin_y>=0 && bin_y<coverage_size)
                        coverage_false[i_coverage][bin_x*coverage_size+bin_y]++;
                  }
               }
            }
            for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
              bin_x=(double)(coverage_size)*(P_last[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
              for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
                bin_y=(double)(coverage_size)*(P_last[j_P]-P_min[j_P])/(P_max[j_P]-P_min[j_P]);
                if(bin_x>=0 && bin_x<coverage_size && bin_y>=0 && bin_y<coverage_size){
                  coverage_keep[i_coverage][bin_x*coverage_size+bin_y]++;
                  max_L_surface[i_coverage][bin_x*coverage_size+bin_y]  =MAX(max_L_surface[i_coverage][bin_x*coverage_size+bin_y],ln_likelihood_last);
                  mean_L_surface[i_coverage][bin_x*coverage_size+bin_y]+=ln_likelihood_last;
                }
              }
            }

            // Build histograms
            for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
              bin_x=(double)(coverage_size)*(P_last[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
              if(bin_x>=0 && bin_x<coverage_size){
                P_histogram[i_P][bin_x]++;
                mean_L_histogram[i_P][bin_x]+=ln_likelihood_last;
                max_L_histogram[i_P][bin_x]  =MAX(max_L_histogram[i_P][bin_x],ln_likelihood_last);
              }
            }
            for(i_DS=0;i_DS<n_DS;i_DS++){
              for(i_M=0;i_M<n_M[i_DS];i_M++){
                bin_x=(double)(coverage_size)*(M_last[i_DS][i_M]-M_min[i_DS][i_M])/(M_max[i_DS][i_M]-M_min[i_DS][i_M]);
                if(bin_x>=0 && bin_x<coverage_size)
                  M_histogram[i_DS][i_M][bin_x]++;
              }
            }
      
            // Compute parameter standard deviations
            for(i_P=0;i_P<n_P;i_P++){
              dP_avg[i_P]+=pow(P_last[i_P]-P_avg[i_P],2.);
            }
            if(!flag_no_map_write){
              for(i_DS=0;i_DS<n_DS;i_DS++){
                for(i_M=0;i_M<n_M[i_DS];i_M++)
                  dM_avg[i_DS][i_M]+=pow(M_last[i_DS][i_M]-M_avg[i_DS][i_M],2.);
              }
            }
          }
        }
      } // i_iteration
      fclose(fp_chain);

      // Finish mean likelihood histograms & surfaces
      for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
        for(j_coverage=0;j_coverage<coverage_size;j_coverage++)
           mean_L_histogram[i_P][j_coverage]/=(double)P_histogram[i_P][j_coverage];
        for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
          for(j_coverage=0;j_coverage<coverage_size*coverage_size;j_coverage++)
            mean_L_surface[i_coverage][j_coverage]/=(double)coverage_keep[i_coverage][j_coverage];
        }
      }

      // Finish standard deviations
      for(i_P=0;i_P<n_P;i_P++)
        dP_avg[i_P]=sqrt(dP_avg[i_P]/(double)n_used);
      if(!flag_no_map_write){
        for(i_DS=0;i_DS<n_DS;i_DS++){
          for(i_M=0;i_M<n_M[i_DS];i_M++)
            dM_avg[i_DS][i_M]=sqrt(dM_avg[i_DS][i_M]/(double)n_used);
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);    

      // Finish covariance matrix
      for(i_P=0;i_P<n_P;i_P++)
        P_i_bar[i_P]=P_i_bar_accum[i_P]/(double)n_used;
      for(i_P=0,k_P=0;i_P<n_P;i_P++){
        for(j_P=0;j_P<n_P;j_P++,k_P++){
          P_ij_bar[k_P] =P_ij_bar_accum[k_P]/(double)n_used;
          V_compute[k_P]=P_ij_bar[k_P]-P_i_bar[i_P]*P_i_bar[j_P];
        }
      }

      // Write covariance matrix
      if(SID.I_am_Master){
        fp_covariance=fopen(filename_covariance,"wb");
        fwrite(&n_P,     sizeof(int),   1,      fp_covariance);
        fwrite(V_compute,sizeof(double),n_P*n_P,fp_covariance);
        fclose(fp_covariance);
      }

      // Compute mapped dataset confidence intervals from histograms
      SID_log("Compute confidence intervals...",SID_LOG_OPEN);
      if(!flag_no_map_write){
        for(i_DS=0;i_DS<n_DS;i_DS++){
          for(i_M=0;i_M<n_M[i_DS];i_M++){
            merge_sort(M_histogram[i_DS][i_M],(size_t)coverage_size,&histogram_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
            accumulator       =M_histogram[i_DS][i_M][histogram_index[coverage_size-1]];
            M_best_index      =histogram_index[coverage_size-1];
            M_lo_68_index     =histogram_index[coverage_size-1];
            M_hi_68_index     =histogram_index[coverage_size-1];
            for(j_P=coverage_size-2;j_P>=0 && ((double)accumulator/(double)n_used)<0.68;j_P--){
              if(histogram_index[j_P]<M_lo_68_index) M_lo_68_index=histogram_index[j_P];
              if(histogram_index[j_P]>M_hi_68_index) M_hi_68_index=histogram_index[j_P];
              accumulator+=M_histogram[i_DS][i_M][histogram_index[j_P]];
            }
            M_lo_95_index=M_lo_68_index;
            M_hi_95_index=M_hi_68_index;
            for(;j_P>=0 && ((double)accumulator/(double)n_used)<0.95;j_P--){
              if(histogram_index[j_P]<M_lo_95_index) M_lo_95_index=histogram_index[j_P];
              if(histogram_index[j_P]>M_hi_95_index) M_hi_95_index=histogram_index[j_P];
              accumulator+=M_histogram[i_DS][i_M][histogram_index[j_P]];
            }
            M_best[i_DS][i_M] =M_min[i_DS][i_M]+(M_max[i_DS][i_M]-M_min[i_DS][i_M])*((double)(M_best_index)+0.5)/(double)coverage_size;
            M_lo_68[i_DS][i_M]=M_min[i_DS][i_M]+(M_max[i_DS][i_M]-M_min[i_DS][i_M])*((double)(M_lo_68_index)+0.5)/(double)coverage_size;
            M_hi_68[i_DS][i_M]=M_min[i_DS][i_M]+(M_max[i_DS][i_M]-M_min[i_DS][i_M])*((double)(M_hi_68_index)+0.5)/(double)coverage_size;
            M_lo_95[i_DS][i_M]=M_min[i_DS][i_M]+(M_max[i_DS][i_M]-M_min[i_DS][i_M])*((double)(M_lo_95_index)+0.5)/(double)coverage_size;
            M_hi_95[i_DS][i_M]=M_min[i_DS][i_M]+(M_max[i_DS][i_M]-M_min[i_DS][i_M])*((double)(M_hi_95_index)+0.5)/(double)coverage_size;
            SID_free(SID_FARG histogram_index);
          }
        }
      }

      // Compute parameter confidence intervals from histograms
      for(i_P=0;i_P<n_P;i_P++){
        merge_sort(P_histogram[i_P],(size_t)coverage_size,&histogram_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
        accumulator       =P_histogram[i_P][histogram_index[coverage_size-1]];
        P_best_index      =histogram_index[coverage_size-1];
        P_lo_68_index     =histogram_index[coverage_size-1];
        P_hi_68_index     =histogram_index[coverage_size-1];
        for(j_P=1,k_P=0;j_P<coverage_size;j_P++)
           if(max_L_histogram[i_P][j_P]>max_L_histogram[i_P][k_P]) k_P=j_P;
        P_peak[i_P]=P_min[i_P]+(P_max[i_P]-P_min[i_P])*((double)(k_P)+0.5)/(double)coverage_size;
        for(j_P=coverage_size-2;j_P>=0 && ((double)accumulator/(double)n_used)<0.68;j_P--){
          if(histogram_index[j_P]<P_lo_68_index) P_lo_68_index=histogram_index[j_P];
          if(histogram_index[j_P]>P_hi_68_index) P_hi_68_index=histogram_index[j_P];
          accumulator+=P_histogram[i_P][histogram_index[j_P]];
        }
        P_lo_95_index=P_lo_68_index;
        P_hi_95_index=P_hi_68_index;
        for(;j_P>=0 && ((double)accumulator/(double)n_used)<0.95;j_P--){
          if(histogram_index[j_P]<P_lo_95_index) P_lo_95_index=histogram_index[j_P];
          if(histogram_index[j_P]>P_hi_95_index) P_hi_95_index=histogram_index[j_P];
          accumulator+=P_histogram[i_P][histogram_index[j_P]];
        }
        P_best[i_P] =P_min[i_P]+(P_max[i_P]-P_min[i_P])*((double)(P_best_index)+0.5)/(double)coverage_size;
        P_lo_68[i_P]=P_min[i_P]+(P_max[i_P]-P_min[i_P])*((double)(P_lo_68_index)+0.5)/(double)coverage_size;
        P_hi_68[i_P]=P_min[i_P]+(P_max[i_P]-P_min[i_P])*((double)(P_hi_68_index)+0.5)/(double)coverage_size;
        P_lo_95[i_P]=P_min[i_P]+(P_max[i_P]-P_min[i_P])*((double)(P_lo_95_index)+0.5)/(double)coverage_size;
        P_hi_95[i_P]=P_min[i_P]+(P_max[i_P]-P_min[i_P])*((double)(P_hi_95_index)+0.5)/(double)coverage_size;
        SID_free(SID_FARG histogram_index);
      }

      // Compute parameter confidence contours from coverage maps
      for(i_coverage=0;i_coverage<n_coverage;i_coverage++){
        merge_sort(coverage_keep[i_coverage],(size_t)(coverage_size*coverage_size),&coverage_keep_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
        accumulator             =coverage_keep[i_coverage][coverage_keep_index[coverage_size*coverage_size-1]];
        P_contour_68[i_coverage]=coverage_keep[i_coverage][coverage_keep_index[coverage_size*coverage_size-1]];
        for(j_P=(coverage_size*coverage_size)-2;j_P>=0 && ((double)accumulator/(double)n_used)<0.68;j_P--){
          P_contour_68[i_coverage] =coverage_keep[i_coverage][coverage_keep_index[j_P]];
          accumulator             +=coverage_keep[i_coverage][coverage_keep_index[j_P]];
        }
        P_contour_95[i_coverage]=P_contour_68[i_coverage];
        for(;j_P>=0 && ((double)accumulator/(double)n_used)<0.95;j_P--){
          P_contour_95[i_coverage] =coverage_keep[i_coverage][coverage_keep_index[j_P]];
          accumulator             +=coverage_keep[i_coverage][coverage_keep_index[j_P]];
        }
        SID_free(SID_FARG coverage_keep_index);
      }
      SID_log("Done.",SID_LOG_CLOSE);

      // Write coverage maps
      SID_log("Writing coverage maps...",SID_LOG_OPEN);
      fp_coverage=fopen(filename_coverage,"wb");
      fwrite(&n_coverage,   sizeof(int),     1,fp_coverage);
      fwrite(&coverage_size,sizeof(int),     1,fp_coverage);
      fwrite(P_min,         sizeof(double),n_P,fp_coverage);
      fwrite(P_max,         sizeof(double),n_P,fp_coverage);
      for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
        for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
          fwrite(&P_contour_68[i_coverage],sizeof(double),1,fp_coverage);
          fwrite(&P_contour_95[i_coverage],sizeof(double),1,fp_coverage);

          // We need to covert to uint64_t so that allresults_MCMC.py works on 32-bit systems
          for(j_coverage=0;j_coverage<coverage_size*coverage_size;j_coverage++)
            temp_out[j_coverage]=coverage_true[i_coverage][j_coverage];
          fwrite(temp_out,sizeof(uint64_t),coverage_size*coverage_size,fp_coverage);
          for(j_coverage=0;j_coverage<coverage_size*coverage_size;j_coverage++)
            temp_out[j_coverage]=coverage_false[i_coverage][j_coverage];
          fwrite(temp_out,sizeof(uint64_t),coverage_size*coverage_size,fp_coverage);
          for(j_coverage=0;j_coverage<coverage_size*coverage_size;j_coverage++)
            temp_out[j_coverage]=coverage_keep[i_coverage][j_coverage];
          fwrite(temp_out,sizeof(uint64_t),coverage_size*coverage_size,fp_coverage);

          fwrite(max_L_surface[i_coverage], sizeof(double),coverage_size*coverage_size,fp_coverage);
          fwrite(mean_L_surface[i_coverage],sizeof(double),coverage_size*coverage_size,fp_coverage);
        }
      }
      fclose(fp_coverage);
      SID_log("Done.",SID_LOG_CLOSE);

      // Write histograms
      SID_log("Writing histograms...",SID_LOG_OPEN);
      fp_histograms=fopen(filename_histograms,"wb");
      for(i_P=0;i_P<n_P;i_P++){
        fwrite(&P_best[i_P],    sizeof(double),1,            fp_histograms);
        fwrite(&P_lo_68[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_hi_68[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_lo_95[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_hi_95[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_peak[i_P],    sizeof(double),1,            fp_histograms);
        // We need to covert to uint64_t so that allresults_MCMC.py works on 32-bit systems
        for(j_coverage=0;j_coverage<coverage_size;j_coverage++)
          temp_out[j_coverage]=P_histogram[i_P][j_coverage];
        fwrite(temp_out,             sizeof(uint64_t),coverage_size,fp_histograms);
        fwrite(mean_L_histogram[i_P],sizeof(double),  coverage_size,fp_histograms);
        fwrite(max_L_histogram[i_P], sizeof(double),  coverage_size,fp_histograms);
      }
      if(!flag_no_map_write){
        for(i_DS=0;i_DS<n_DS;i_DS++){
          for(i_M=0;i_M<n_M[i_DS];i_M++){
            fwrite(&M_best[i_DS][i_M],    sizeof(double),1,            fp_histograms);
            fwrite(&M_lo_68[i_DS][i_M],   sizeof(double),1,            fp_histograms);
            fwrite(&M_hi_68[i_DS][i_M],   sizeof(double),1,            fp_histograms);
            fwrite(&M_lo_95[i_DS][i_M],   sizeof(double),1,            fp_histograms);
            fwrite(&M_hi_95[i_DS][i_M],   sizeof(double),1,            fp_histograms);
            // We need to covert to uint64_t so that allresults_MCMC.py works on 32-bit systems
            for(j_coverage=0;j_coverage<coverage_size;j_coverage++)
              temp_out[j_coverage]=M_histogram[i_DS][i_M][j_coverage];
            fwrite(temp_out,sizeof(uint64_t),coverage_size,fp_histograms);
          }
        }
      }
      fclose(fp_histograms);
      SID_log("Done.",SID_LOG_CLOSE);

      // Compute best fit models and their likelihoods
      if(MCMC->map_P_to_M!=NULL){
        MCMC->map_P_to_M(P_best,MCMC,M_best_parameters);
        MCMC->compute_MCMC_ln_likelihood(MCMC,M_best_parameters,P_best,&ln_likelihood_best);
        MCMC->map_P_to_M(P_peak,MCMC,M_peak_parameters);
        MCMC->compute_MCMC_ln_likelihood(MCMC,M_peak_parameters,P_peak,&ln_likelihood_peak);
      }

      // Write fit results
      SID_log("Writing final statistics...",SID_LOG_OPEN);
      sprintf(filename_results,"%s/fit_for_parameters.dat",filename_results_dir);
      fp_results=fopen(filename_results,"w");
      fprintf(fp_results,"# MCMC parameter fit results to %s\n",MCMC->problem_name);
      if(MCMC->map_P_to_M!=NULL){
        fprintf(fp_results,"#   ln(likelihood)+constant=%le (marginalized PDF)\n",  ln_likelihood_best);
        fprintf(fp_results,"#   ln(likelihood)+constant=%le (maximum likelihood)\n",ln_likelihood_peak);
      }
      fprintf(fp_results,"#   n_samples_used=%d\n",n_used);
      fprintf(fp_results,"#   n_parameters  =%d\n",n_P);
      fprintf(fp_results,"#   n_data_sets   =%d\n",n_DS);
      fprintf(fp_results,"#   used data sets:\n");
      current_DS=MCMC->DS;
      while(current_DS!=NULL){
        next_DS=current_DS->next;
        fprintf(fp_results,"#     %s\n",current_DS->name);
        current_DS=next_DS;
        i_DS++;
      }
      sprintf(column_txt,"Column:");
      i_column=1;
      fprintf(fp_results,"# %s (%02d) Parameter name\n",               column_txt,i_column++);
      sprintf(column_txt,"       ");
      for(i_array=0;i_array<MCMC->n_arrays;i_array++)
        fprintf(fp_results,"# %s (%02d) %s\n",                         column_txt,i_column++,MCMC->array_name[i_array]);
      fprintf(fp_results,"# %s (%02d) Initial value\n",                column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Average\n",                      column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Standard deviation\n",           column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Maximum likelihood value\n",     column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Lower limit (68%% confidence)\n",column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Upper limit (68%% confidence)\n",column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Lower limit (95%% confidence)\n",column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Upper limit (95%% confidence)\n",column_txt,i_column++);
      fprintf(fp_results,"# %s (%02d) Best fit value\n",               column_txt,i_column++);
      for(i_P=0;i_P<n_P;i_P++){
        strcpy(array_name_temp,MCMC->P_names[i_P]);
        search_and_replace(array_name_temp," ","~");
        fprintf(fp_results,MCMC->P_name_format,array_name_temp);
        for(i_array=0;i_array<MCMC->n_arrays;i_array++)
          fprintf(fp_results," %13.6le",MCMC->array[i_array][i_P]);
        fprintf(fp_results," %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le\n",MCMC->P_init[i_P],P_avg[i_P],dP_avg[i_P],P_best[i_P],P_lo_68[i_P],P_hi_68[i_P],P_lo_95[i_P],P_hi_95[i_P],P_peak[i_P]);
      }
      fclose(fp_results);

      // Write the fits for each dataset
      i_DS=0;
      current_DS=MCMC->DS;
      while(current_DS!=NULL){
        next_DS  =current_DS->next;
        M_target =current_DS->M_target;
        dM_target=current_DS->dM_target;
        if(current_DS->n_D<=1){
          if(n_DS>1)
            sprintf(filename_results,"%s/fit_for_dataset_%05d.dat",filename_results_dir,i_DS);
          else
            sprintf(filename_results,"%s/fit_for_dataset.dat",filename_results_dir);
          fp_results=fopen(filename_results,"w");
          fprintf(fp_results,"# MCMC data set fit results for %s\n",current_DS->name);
          fprintf(fp_results,"#   n_samples_used=%d\n",n_used);
          fprintf(fp_results,"#   dataset_DoF   =%d\n",n_M[i_DS]);
          sprintf(column_txt,"Column:");
          for(i_array=0,i_column=1;i_array<MCMC->n_arrays;i_array++){
            fprintf(fp_results,"# %s (%02d) %s\n",                         column_txt,i_column++,current_DS->array_name[i_array]);
            sprintf(column_txt,"       ");
          }
          sprintf(column_txt,"       ");
          fprintf(fp_results,"# %s (%02d) Dataset\n",                      column_txt,i_column++);
          fprintf(fp_results,"# %s (%02d) Dataset uncertainty\n",          column_txt,i_column++);
          if(!flag_no_map_write){
            fprintf(fp_results,"# %s (%02d) Average\n",                      column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Standard deviation\n",           column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Marginalized likelihood value\n",column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Lower limit (68%% confidence)\n",column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Upper limit (68%% confidence)\n",column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Lower limit (95%% confidence)\n",column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Upper limit (95%% confidence)\n",column_txt,i_column++);
          }
          if(MCMC->map_P_to_M!=NULL){
            fprintf(fp_results,"# %s (%02d) Model (Marginalized likelihood parameters)\n",column_txt,i_column++);
            fprintf(fp_results,"# %s (%02d) Model (Maximum      likelihood parameters)\n",column_txt,i_column++);
          }
          for(i_P=0;i_P<n_M[i_DS];i_P++){
            for(i_array=0;i_array<n_M_arrays[i_DS];i_array++)
              fprintf(fp_results,"%13.6le ",M_arrays[i_DS][i_array][i_P]);
            fprintf(fp_results,"%13.6le %13.6le",M_target[i_P],dM_target[i_P]);
            if(!flag_no_map_write)
              fprintf(fp_results," %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le",
                      M_avg[i_DS][i_P],
                      dM_avg[i_DS][i_P],
                      M_best[i_DS][i_P],
                      M_lo_68[i_DS][i_P],
                      M_hi_68[i_DS][i_P],
                      M_lo_95[i_DS][i_P],
                      M_hi_95[i_DS][i_P]);
            if(MCMC->map_P_to_M!=NULL)
              fprintf(fp_results,"  %13.6le %13.6le",
                      M_best_parameters[i_DS][i_P],
                      M_peak_parameters[i_DS][i_P]);
            fprintf(fp_results,"\n");
          }
          fclose(fp_results);
        }
        else{
          if(n_DS>1)
            sprintf(filename_results,"%s/fit_for_dataset_%05d.fits",filename_results_dir,i_DS);
          else
            sprintf(filename_results,"%s/fit_for_dataset.fits",filename_results_dir);
          write_image_FITS(  M_target,SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"DATASET");
          append_image_FITS(dM_target,SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"UNCERTAINTY");
          if(!flag_no_map_write){
            append_image_FITS(M_best[i_DS], SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MAP_MAXL");
            append_image_FITS(M_lo_68[i_DS],SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MAP_68PC_LO");
            append_image_FITS(M_hi_68[i_DS],SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MAP_68PC_HI");
            append_image_FITS(M_lo_95[i_DS],SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MAP_95PC_LO");
            append_image_FITS(M_hi_95[i_DS],SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MAP_95PC_HI");
          }
          if(MCMC->map_P_to_M!=NULL){
            // Write the fit from the maximum likelihood of the marginalized pdf
            append_image_FITS(M_best_parameters[i_DS],SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MODEL_MAXL");
            for(i_P=0;i_P<n_M[i_DS];i_P++) residual[i_P]=M_target[i_P]-M_best_parameters[i_DS][i_P];
            append_image_FITS(residual,SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"RESIDUAL_MAXL");
            // Write the best fit
            append_image_FITS(M_peak_parameters[i_DS],SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"MODEL_BEST");
            for(i_P=0;i_P<n_M[i_DS];i_P++) residual[i_P]=M_target[i_P]-M_peak_parameters[i_DS][i_P];
            append_image_FITS(residual,SID_DOUBLE,current_DS->n_D,current_DS->D,filename_results,"RESIDUAL_BEST");
          }
        }
        current_DS=next_DS;
        i_DS++;
      }
      SID_log("Done.",SID_LOG_CLOSE);
  } // if this is my chain

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  SID_free(SID_FARG P_i_bar_accum);
  SID_free(SID_FARG P_ij_bar_accum);
  SID_free(SID_FARG P_i_bar);
  SID_free(SID_FARG P_ij_bar);
  SID_free(SID_FARG V_compute);
  for(i_DS=0;i_DS<n_DS;i_DS++){
    SID_free(SID_FARG M_avg[i_DS]);
    SID_free(SID_FARG dM_avg[i_DS]);
    SID_free(SID_FARG M_best[i_DS]);
    SID_free(SID_FARG M_best_parameters[i_DS]);
    SID_free(SID_FARG M_peak_parameters[i_DS]);
    SID_free(SID_FARG M_min[i_DS]);
    SID_free(SID_FARG M_max[i_DS]);
    SID_free(SID_FARG M_lo_68[i_DS]);
    SID_free(SID_FARG M_hi_68[i_DS]);
    SID_free(SID_FARG M_lo_95[i_DS]);
    SID_free(SID_FARG M_hi_95[i_DS]);
    for(j_M=0;j_M<n_M[i_DS];j_M++)
      SID_free(SID_FARG M_histogram[i_DS][j_M]);
    SID_free(SID_FARG M_histogram[i_DS]);
  }
  SID_free(SID_FARG M_avg);
  SID_free(SID_FARG dM_avg);
  SID_free(SID_FARG M_best);
  SID_free(SID_FARG residual);
  SID_free(SID_FARG M_best_parameters);
  SID_free(SID_FARG M_peak_parameters);
  SID_free(SID_FARG M_min);
  SID_free(SID_FARG M_max);
  SID_free(SID_FARG M_lo_68);
  SID_free(SID_FARG M_hi_68);
  SID_free(SID_FARG M_lo_95);
  SID_free(SID_FARG M_hi_95);
  SID_free(SID_FARG M_histogram);
  SID_free(SID_FARG P_contour_68);
  SID_free(SID_FARG P_contour_95);
  for(i_P=0;i_P<n_P;i_P++){
    SID_free(SID_FARG P_histogram[i_P]);
    SID_free(SID_FARG mean_L_histogram[i_P]);
    SID_free(SID_FARG max_L_histogram[i_P]);
  }
  SID_free(SID_FARG P_histogram);
  SID_free(SID_FARG n_M_arrays);
  SID_free(SID_FARG M_arrays);
  SID_free(SID_FARG M_avg);
  SID_free(SID_FARG dM_avg);
  for(i_P=0;i_P<n_coverage;i_P++){
    SID_free(SID_FARG coverage_true[i_P]);
    SID_free(SID_FARG coverage_false[i_P]);
    SID_free(SID_FARG coverage_keep[i_P]);
    SID_free(SID_FARG max_L_surface[i_P]);
    SID_free(SID_FARG mean_L_surface[i_P]);
  }
  SID_free(SID_FARG coverage_true);
  SID_free(SID_FARG coverage_false);
  SID_free(SID_FARG coverage_keep);
  SID_free(SID_FARG max_L_surface);
  SID_free(SID_FARG mean_L_surface);
  SID_free(SID_FARG temp_out);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}
