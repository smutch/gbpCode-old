#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void compute_MCMC(MCMC_info *MCMC){
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
  int       n_avg_test;
  int       n_avg_covariance_test;
  int       flag_autocor_on_test;
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
  double    ln_Pr_last;
  double    ln_Pr_new;
  double    ln_Pr_min;
  double    ln_Pr_avg;
  double    ln_Pr_max;
  double  **M_best;
  double  **M_min;
  double  **M_max;
  double  **M_new;
  double  **M_lo_68;
  double  **M_hi_68;
  double  **M_lo_95;
  double  **M_hi_95;
  size_t ***M_histogram;
  double  **M_last;
  double  **M_avg;
  double  **dM_avg;
  double   *slopes;
  double   *dP_sub;
  double   *ln_Pr_chain;
  double   *ln_Pr_proposals;
  double   *drift;
  double  **auto_cor;
  double  **P_chain;
  double  **P_proposals;
  size_t **P_histogram;
  size_t **coverage_true;
  size_t **coverage_false;
  size_t **coverage_keep;
  double    L_x,L_y,L_z;
  int       n_x,n_y,n_z;
  int       n_C_x;
  double    ln_likelihood_new,ln_likelihood_last;
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
  RNG_info  RNG;
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
  int       flag_report_props;
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
  int             n_avg_covariance;
  int             n_iterations_file_total;
  int             n_iterations_file_burn;
  int             flag_restart=FALSE;

  SID_log("Performing MCMC...",SID_LOG_OPEN|SID_LOG_TIMER);

  SID_log("Initializing MCMC...",SID_LOG_OPEN);

  n_P                   =MCMC->n_P;
  n_DS                  =MCMC->n_DS;
  temperature           =MCMC->temperature;
  n_avg                 =MCMC->n_avg;
  n_avg_covariance      =MCMC->n_avg_covariance;
  n_iterations          =MCMC->n_iterations;
  n_iterations_burn     =MCMC->n_iterations_burn;
  n_iterations_integrate=n_iterations-n_iterations_burn;
  n_thin                =MCMC->n_thin;
  coverage_size         =MCMC->coverage_size;
  flag_autocor_on       =MCMC->flag_autocor_on;
  V                     =MCMC->V;
  m                     =MCMC->m;
  b                     =gsl_vector_calloc(n_P);

  flag_report_props     =check_mode_for_flag(MCMC->mode,MCMC_REPORT_PROPS);
  
  strcpy(filename_output_dir,MCMC->filename_output_dir);

  constant=MCMC->temperature*2.4/sqrt((double)MCMC->n_M_total);

  SID_log("temperature            = %lf",SID_LOG_COMMENT,temperature);
  SID_log("n_avg                  = %d",SID_LOG_COMMENT,n_avg);
  SID_log("n_avg_covariance       = %d",SID_LOG_COMMENT,n_avg_covariance);
  SID_log("n_thin                 = %d",SID_LOG_COMMENT,n_thin);
  SID_log("coverage_size          = %d",SID_LOG_COMMENT,coverage_size);
  SID_log("n_iterations           = %d",SID_LOG_COMMENT,n_iterations);
  SID_log("n_iterations_burn      = %d",SID_LOG_COMMENT,n_iterations_burn);
  SID_log("n_iterations_integrate = %d",SID_LOG_COMMENT,n_iterations_integrate);

  // Initialize random number generator
  init_seed_from_clock(&seed);
  SID_log("random seed            = %d",SID_LOG_COMMENT,seed);
  init_RNG(&seed,&RNG,RNG_DEFAULT);

  if(V!=NULL)
    SID_log("Covariance matrix is initialized.",SID_LOG_COMMENT);
  else
    SID_log("Covariance matrix is NOT initialized.",SID_LOG_COMMENT);
  if(flag_autocor_on)
    SID_log("Auto-correlation  is on.",SID_LOG_COMMENT);
  else  
    SID_log("Auto-correlation  is off.",SID_LOG_COMMENT);

  // Initialize arrays used for computing the covariance matrix
  P_i_bar_accum =(double *)malloc(sizeof(double)*n_P);
  P_ij_bar_accum=(double *)malloc(sizeof(double)*n_P*n_P);
  P_i_bar       =(double *)malloc(sizeof(double)*n_P);
  P_ij_bar      =(double *)malloc(sizeof(double)*n_P*n_P);
  V_compute     =(double *)malloc(sizeof(double)*n_P*n_P);
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
  if(!check_mode_for_flag(MCMC->mode,MCMC_NO_CVM_UPDATE))
    add_covariance_to_MCMC(MCMC,V_compute);

  // Initialize coverage arrays
  for(i_P=0,n_coverage=0;i_P<n_P;i_P++){
    for(j_P=i_P+1;j_P<n_P;j_P++)
      n_coverage++;
  }
  coverage_true =(size_t **)SID_malloc(sizeof(size_t *)*n_coverage);
  coverage_false=(size_t **)SID_malloc(sizeof(size_t *)*n_coverage);
  coverage_keep =(size_t **)SID_malloc(sizeof(size_t *)*n_coverage);
  P_contour_68  =(double *)SID_malloc(sizeof(double)*n_coverage);
  P_contour_95  =(double *)SID_malloc(sizeof(double)*n_coverage);
  for(i_P=0;i_P<n_coverage;i_P++){
    coverage_true[i_P] =(size_t *)SID_malloc(sizeof(size_t)*coverage_size*coverage_size);
    coverage_false[i_P]=(size_t *)SID_malloc(sizeof(size_t)*coverage_size*coverage_size);
    coverage_keep[i_P] =(size_t *)SID_malloc(sizeof(size_t)*coverage_size*coverage_size);
  }
  P_histogram=(size_t **)SID_malloc(sizeof(size_t *)*n_P);
  for(i_P=0;i_P<n_P;i_P++)
    P_histogram[i_P]   =(size_t *)SID_malloc(sizeof(size_t)*coverage_size);

  // Initialize chain arrays
  n_M        =(int      *)SID_malloc(sizeof(int      )*n_DS);
  M_best     =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_min      =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_max      =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_lo_68    =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_hi_68    =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_lo_95    =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_hi_95    =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_new      =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_last     =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_avg      =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  dM_avg     =(double  **)SID_malloc(sizeof(double  *)*n_DS);
  M_histogram=(size_t ***)SID_malloc(sizeof(size_t **)*n_DS);
  current_DS=MCMC->DS;
  i_DS      =0;
  while(current_DS!=NULL){
    next_DS              =current_DS->next;
    n_M[i_DS]            =current_DS->n_M;
    M_best[i_DS]     =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_min[i_DS]      =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_max[i_DS]      =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_lo_68[i_DS]    =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_hi_68[i_DS]    =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_lo_95[i_DS]    =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_hi_95[i_DS]    =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_best[i_DS]     =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_new[i_DS]      =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_last[i_DS]     =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_avg[i_DS]      =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    dM_avg[i_DS]     =(double *)SID_malloc(sizeof(double)*n_M[i_DS]);
    M_histogram[i_DS]=(size_t **)SID_malloc(sizeof(size_t *)*n_M[i_DS]);
    for(j_M=0;j_M<n_M[i_DS];j_M++)
      M_histogram[i_DS][j_M]=(size_t *)SID_malloc(sizeof(size_t)*coverage_size);
    current_DS=next_DS;
    i_DS++;
  }
  P_best         =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_min          =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_max          =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_new          =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_last         =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_chain        =(double **)SID_malloc(sizeof(double *)*n_P);
  P_proposals    =(double **)SID_malloc(sizeof(double *)*n_P);
  P_avg          =(double  *)SID_malloc(sizeof(double)  *n_P);
  dP_avg         =(double  *)SID_malloc(sizeof(double)  *n_P);
  slopes         =(double  *)SID_malloc(sizeof(double *)*n_P);
  dP_sub         =(double  *)SID_malloc(sizeof(double *)*n_P);
  drift          =(double  *)SID_malloc(sizeof(double *)*n_P);
  P_lo_68        =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_hi_68        =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_lo_95        =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_hi_95        =(double  *)SID_malloc(sizeof(double)  *n_P);
  ln_Pr_chain    =(double  *)SID_malloc(sizeof(double)*n_avg);
  ln_Pr_proposals=(double  *)SID_malloc(sizeof(double)*n_avg);
  for(i_P=0;i_P<n_P;i_P++)
    P_last[i_P]=MCMC->P_init[i_P];
  if(flag_autocor_on)
    auto_cor=(double **)SID_malloc(sizeof(double *)*n_P);
  else
    auto_cor=NULL;
  for(i_P=0;i_P<n_P;i_P++){
    P_chain[i_P]    =(double *)SID_malloc(sizeof(double)*n_avg);
    P_proposals[i_P]=(double *)SID_malloc(sizeof(double)*n_avg);
    if(flag_autocor_on)
      auto_cor[i_P]=(double *)SID_malloc(sizeof(double)*n_avg);
  }
  for(i_P=0;i_P<n_P;i_P++)
    drift[i_P]=0.;

  // Create arrays of pointers for the associate P&M arrays for quicker/easier look-up
  n_M_arrays=(int      *)SID_malloc(sizeof(int      )*n_DS);
  M_arrays  =(double ***)SID_malloc(sizeof(double **)*n_DS);
  current_DS=MCMC->DS;
  i_DS      =0;
  while(current_DS!=NULL){
    next_DS         =current_DS->next;
    n_M_arrays[i_DS]=current_DS->n_arrays;
    M_arrays[i_DS]  =current_DS->array;
    current_DS      =next_DS;
    i_DS++;
  }

  // Set the local chain number
  if(check_mode_for_flag(MCMC->mode,MCMC_MODE_PARALLEL))
    my_chain=SID.My_rank;
  else
    my_chain=MASTER_RANK;

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
  
  SID_log("Done.",SID_LOG_CLOSE);

  // Perform integration
  if(MCMC->flag_integrate_on){
    
    SID_log("Initializing integration...",SID_LOG_OPEN);
    
    // Make sure the needed directories exist
    mkdir(filename_output_dir, 02755);
    mkdir(filename_chain_dir,  02755);
    mkdir(filename_results_dir,02755);
    mkdir(filename_plots_dir,  02755);

    // Read/Write Header file
    if(SID.I_am_Master){
      // If the header file already exists, check that it is consistant with the current run
      if((fp_run=fopen(filename_run,"rb"))!=NULL){
        flag_restart=TRUE;
        SID_log("Checking the consistancy of this run with the previous run...",SID_LOG_OPEN);
        fp_run=fopen(filename_run,"rb");
        fread(problem_name_test,sizeof(char),MCMC_NAME_SIZE,fp_run);
        if(strcmp(problem_name_test,MCMC->problem_name))
          SID_trap_error("Problem names are inconsistant (i.e. {%s}!={%s}).",ERROR_LOGIC,MCMC->problem_name,problem_name_test);
        fread(&n_avg_test,sizeof(int),1,fp_run);
        if(n_avg_test!=n_avg)
          SID_trap_error("Integration averaging intervals are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_avg,n_avg_test);
        fread(&n_avg_covariance_test,sizeof(int),1,fp_run);
        if(n_avg_covariance_test!=n_avg_covariance)
          SID_trap_error("Covariance averaging intervals are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_avg_covariance,n_avg_covariance_test);
        fread(&flag_autocor_on_test,sizeof(int),1,fp_run);
        if(flag_autocor_on_test!=flag_autocor_on)
          SID_trap_error("Autocorrelation flags are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,flag_autocor_on,flag_autocor_on_test);
        fread(&n_P_test,sizeof(int),1,fp_run);
        if(n_P_test!=n_P)
          SID_trap_error("The number of paramaters is inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_P,n_P_test);
        for(i_P=0;i_P<n_P;i_P++){
          fread(P_name_test, sizeof(char),  MCMC_NAME_SIZE,fp_run);
          if(strcmp(P_name_test,MCMC->P_names[i_P]))
            SID_trap_error("Parameter #%d's names are inconsistant (i.e. {%s}!={%s}).",ERROR_LOGIC,i_P,MCMC->P_names[i_P],P_name_test);
          fread(&P_init_test,sizeof(double),1,             fp_run);
          if(P_init_test!=MCMC->P_init[i_P])
            SID_trap_error("Parameter #%d's initial values are inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_P,MCMC->P_init[i_P],P_init_test);
        }
        fread(&n_arrays_test,sizeof(int),1,fp_run);
        if(n_arrays_test!=MCMC->n_arrays)
          SID_trap_error("Numbers of project arrays are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,MCMC->n_arrays,n_arrays_test);
        for(i_array=0;i_array<MCMC->n_arrays;i_array++){
          fread(&array_name_test,sizeof(char),MCMC_NAME_SIZE,fp_run);
          if(strcmp(array_name_test,MCMC->array_name[i_array]))
            SID_trap_error("Project array names are inconsisitant (i.e. {%s}!={%s}).",ERROR_LOGIC,MCMC->array_name[i_array],array_name_test);
        }
        fread(&n_DS_test,sizeof(int),   1,   fp_run);
        if(n_DS_test!=n_DS)
          SID_trap_error("The number of datasets is inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_DS,n_DS_test);
        for(i_DS=0;i_DS<n_DS;i_DS++){
          fread(&n_M_test,sizeof(int),1,fp_run);
          if(n_M_test!=n_M[i_DS])
            SID_trap_error("The sizes of dataset #%d are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,i_DS,n_M_test,n_M[i_DS]);
        }
        current_DS=MCMC->DS;
        i_DS=0;
        while(current_DS!=NULL){
          next_DS=current_DS->next;
          fread(name_test,     sizeof(char),  MCMC_NAME_SIZE, fp_run);
          if(strcmp(name_test,current_DS->name))
            SID_trap_error("Dataset #%d's names are inconsistant (i.e. {%s}!={%s}).",ERROR_LOGIC,i_DS,current_DS->name,name_test);
          for(i_M=0;i_M<current_DS->n_M;i_M++){
            fread(&M_target_test, sizeof(double),1,fp_run);
            if(M_target_test!=current_DS->M_target[i_M])
              SID_trap_error("Dataset #%d, value element #%d is inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_DS,i_M,current_DS->M_target[i_M],M_target_test);
          }
          for(i_M=0;i_M<current_DS->n_M;i_M++){
            fread(&dM_target_test,sizeof(double),1,fp_run);
            if(dM_target_test!=current_DS->dM_target[i_M])
              SID_trap_error("Dataset #%d, uncertainty element #%d is inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_DS,i_M,current_DS->M_target[i_M],M_target_test);
          }
          fread(&n_arrays_test,sizeof(int),1,fp_run);
          if(n_arrays_test!=current_DS->n_arrays)
            SID_trap_error("The number arrays in dataset #%d are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,i_DS,n_arrays_test,current_DS->n_arrays);
          for(i_array=0;i_array<current_DS->n_arrays;i_array++){
            fread(array_name_test,sizeof(char),MCMC_NAME_SIZE,fp_run);
            if(strcmp(array_name_test,current_DS->array_name[i_array]))
              SID_trap_error("Array name #%d for dataset #%d is inconsisitant (i.e. {%s}!={%s}).",ERROR_LOGIC,i_array,i_DS,current_DS->array_name[i_array],array_name_test);
          }
          current_DS=next_DS;
          i_DS++;
        }
        fclose(fp_run);
        SID_log("Done.",SID_LOG_CLOSE);
      }
      // ... else if the header file does not exist, write it
      else{
        SID_log("Write header file...",SID_LOG_OPEN);
        fp_run=fopen(filename_run,"wb");
        // Stuff relating to our MCMC project
        fwrite(MCMC->problem_name,sizeof(char),MCMC_NAME_SIZE,fp_run);
        fwrite(&n_avg,            sizeof(int),   1,   fp_run);
        fwrite(&n_avg_covariance, sizeof(int),   1,   fp_run);
        fwrite(&flag_autocor_on,  sizeof(int),   1,   fp_run);
        fwrite(&n_P,              sizeof(int),   1,   fp_run);
        for(i_P=0;i_P<n_P;i_P++){
          fwrite(MCMC->P_names[i_P],  sizeof(char),  MCMC_NAME_SIZE,fp_run);
          fwrite(&(MCMC->P_init[i_P]),sizeof(double),1,             fp_run);
        }
        fwrite(&(MCMC->n_arrays),sizeof(int),1,fp_run);
        for(i_array=0;i_array<MCMC->n_arrays;i_array++)
          fwrite(MCMC->array_name[i_array],sizeof(char),MCMC_NAME_SIZE,fp_run);
        // Stuff relating to the constraining datasets
        fwrite(&n_DS,sizeof(int),   1,   fp_run);
        fwrite(n_M,  sizeof(int),   n_DS,fp_run);
        current_DS=MCMC->DS;
        while(current_DS!=NULL){
          next_DS=current_DS->next;
          fwrite(current_DS->name,     sizeof(char),  MCMC_NAME_SIZE, fp_run);
          fwrite(current_DS->M_target, sizeof(double),current_DS->n_M,fp_run);
          fwrite(current_DS->dM_target,sizeof(double),current_DS->n_M,fp_run);
          fwrite(&(current_DS->n_arrays),sizeof(int),1,fp_run);
          for(i_array=0;i_array<current_DS->n_arrays;i_array++)
            fwrite(current_DS->array_name[i_array],sizeof(char),MCMC_NAME_SIZE,fp_run);
          current_DS=next_DS;
        }
        fclose(fp_run);
        SID_log("Done.",SID_LOG_CLOSE);
      }
    }
    SID_Bcast(&flag_restart,sizeof(int),MASTER_RANK);

    // Initialize the run parameters
    n_iterations_file_total=0;
    n_iterations_file_burn =0;
    n_covariance           =0;
    i_covariance           =0;

    // If this is NOT a restart, start from scratch ...
    if(!flag_restart){
      SID_log("Setting initial state...",SID_LOG_OPEN);
      fp_chain           =fopen(filename_chain,"wb");
      fp_stats           =fopen(filename_stats,"wb");
      MCMC->map_P_to_M(P_last,MCMC,M_last);
      MCMC->first_map_call=FALSE;
      MCMC->compute_MCMC_ln_likelihood(MCMC,M_last,&ln_likelihood_last);
      fp_chain_iterations=fopen(filename_chain_iterations,"wb");
      fwrite(&n_iterations_file_total,sizeof(int),1,fp_chain_iterations);
      fwrite(&n_iterations_file_burn, sizeof(int),1,fp_chain_iterations);
      fclose(fp_chain_iterations);
      SID_log("Done.",SID_LOG_CLOSE);
    }
    // ... else read the state we left-off from ...
    else if(my_chain==SID.My_rank){
      SID_log("Loading previous state...",SID_LOG_OPEN);
      // ... fetch the number of intervals that have already been computed ...
      SID_log("Reading the existant number of iterations...",SID_LOG_OPEN);
      fp_chain_iterations=fopen(filename_chain_iterations,"rb");
      fread(&n_iterations_file_total,sizeof(int),1,fp_chain_iterations);
      fread(&n_iterations_file_burn, sizeof(int),1,fp_chain_iterations);
      SID_log("# burn  iterations = %d (%d requested)",SID_LOG_COMMENT,n_iterations_file_burn ,n_iterations_burn);
      SID_log("# total iterations = %d (%d requested)",SID_LOG_COMMENT,n_iterations_file_total,n_iterations);
      fclose(fp_chain_iterations);
      SID_log("Done.",SID_LOG_CLOSE);

      // ... read the covariance matrix (if it exists) ...
      if((fp_chain_covariance=fopen(filename_chain_covariance,"rb"))!=NULL){
        SID_log("Reading the existant covariance matrix...",SID_LOG_OPEN);
        fread(&i_covariance,sizeof(int),1,fp_chain_covariance);
        fread(&n_covariance,sizeof(int),1,fp_chain_covariance);
        SID_log("# iterations = %d (averaging interval=%d)",SID_LOG_COMMENT,i_covariance,n_avg_covariance);
        SID_log("# samples    = %d",                        SID_LOG_COMMENT,n_covariance);
        fread(&n_P_test,    sizeof(int),1,fp_chain_covariance);
        if(n_P_test!=n_P)
          SID_trap_error("The number of paramaters is inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_P,n_P_test);        
        fread(V_compute,     sizeof(double),n_P*n_P,fp_chain_covariance);
        fread(P_i_bar_accum, sizeof(double),n_P,    fp_chain_covariance);
        fread(P_ij_bar_accum,sizeof(double),n_P*n_P,fp_chain_covariance);
        fclose(fp_chain_covariance);
        SID_log("Done.",SID_LOG_CLOSE);
        if(!check_mode_for_flag(MCMC->mode,MCMC_NO_CVM_UPDATE))
          add_covariance_to_MCMC(MCMC,V_compute);
      }

      // ... scan to the end of the chain and read the last-used parameter set
      if(n_iterations_file_total>0){
        SID_log("Reading the last-used parameter set...",SID_LOG_OPEN);
        fp_chain=fopen(filename_chain,"rb");
        for(i_iteration=0;i_iteration<n_iterations_file_total;i_iteration++){
          for(i_avg=0;i_avg<n_avg;i_avg++){
            fread(&flag_success,sizeof(char),  1,fp_chain);
            fread(&ln_Pr_last,  sizeof(double),1,fp_chain);
            if(flag_success){
              fread(P_last,sizeof(double),n_P,fp_chain);
              for(i_DS=0;i_DS<n_DS;i_DS++)
                fread(M_last[i_DS],sizeof(double),n_M[i_DS],fp_chain);
            }
            else{            
              fread(P_new,sizeof(double),n_P,fp_chain);
              for(i_DS=0;i_DS<n_DS;i_DS++)
                fread(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
            }
          }        
        }
        fclose(fp_chain);
        SID_log("Done.",SID_LOG_CLOSE);
      }
      else{
        SID_log("Setting initial state...",SID_LOG_OPEN);
        MCMC->map_P_to_M(P_last,MCMC,M_last);
        MCMC->first_map_call=FALSE;
        SID_log("Done.",SID_LOG_CLOSE);
      }
      MCMC->compute_MCMC_ln_likelihood(MCMC,M_last,&ln_likelihood_last);
      
      // ... scan to the end of the chain stats and read the last-generated statistics (particularly the drift)
      SID_log("Reading the last-generated statistics...",SID_LOG_OPEN);
      fp_stats=fopen(filename_stats,"rb");
      for(i_iteration=0;i_iteration<n_iterations_file_total;i_iteration++){
        fread(P_min, sizeof(double),n_P,fp_stats);    // Read proposal stats
        fread(P_avg, sizeof(double),n_P,fp_stats);
        fread(P_max, sizeof(double),n_P,fp_stats);
        fread(dP_avg,sizeof(double),n_P,fp_stats);
        fread(dP_sub,sizeof(double),n_P,fp_stats);
        fread(&ln_Pr_min, sizeof(double),1,fp_stats);
        fread(&ln_Pr_avg, sizeof(double),1,fp_stats);
        fread(&ln_Pr_max, sizeof(double),1,fp_stats);
        if(flag_autocor_on)
          fread(auto_cor,sizeof(double),n_avg-1,fp_stats);
        fread(P_min, sizeof(double),n_P,fp_stats);    // Read chain stats
        fread(P_avg, sizeof(double),n_P,fp_stats);
        fread(P_max, sizeof(double),n_P,fp_stats);
        fread(dP_avg,sizeof(double),n_P,fp_stats);
        fread(dP_sub,sizeof(double),n_P,fp_stats);
        fread(&ln_Pr_min, sizeof(double),1,fp_stats);
        fread(&ln_Pr_avg, sizeof(double),1,fp_stats);
        fread(&ln_Pr_max, sizeof(double),1,fp_stats);
        if(flag_autocor_on)
          fread(auto_cor,sizeof(double),n_avg-1,fp_stats);
        fread(slopes,sizeof(double),n_P,fp_stats);
        fread(drift, sizeof(double),n_P,fp_stats);
      }
      fclose(fp_stats);
      SID_log("Done.",SID_LOG_CLOSE);

      SID_log("Done.",SID_LOG_CLOSE);

      // ... open files ...
      fp_chain=fopen(filename_chain,"ab");
      fp_stats=fopen(filename_stats,"ab");        
    }
    if(flag_report_props && my_chain==SID.My_rank){
      SID_log("Initial parameters:",SID_LOG_ALLRANKS|SID_LOG_OPEN);
      for(i_P=0;i_P<n_P;i_P++)
        SID_log("%s = %le",SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->P_names[i_P],P_last[i_P]);
      SID_log("ln(likelihood)=%le+constant",SID_LOG_ALLRANKS|SID_LOG_COMMENT,ln_likelihood_last);
      SID_log("",SID_LOG_ALLRANKS|SID_LOG_NOPRINT|SID_LOG_CLOSE);
    }
    
    // Remove the existant iterations from the totals we need to perform still
    if(n_iterations_file_total<n_iterations_burn){
      i_phase    =0;
      i_iteration=n_iterations_file_total;
    }
    else if(n_iterations_file_total<n_iterations){
      i_phase    =1;
      i_iteration=n_iterations_file_total-n_iterations_burn;
    }
    else{
      i_phase=2; // ... in other words, we're done; skip the integration.
      i_iteration=n_iterations;
    }
    i_iteration_start=i_iteration;
    SID_log("Done.",SID_LOG_CLOSE);

    // Create the chain in 2 phases: a burn-in phase and an integration phase
    SID_log("Performing integration...",SID_LOG_OPEN|SID_LOG_TIMER);
    for(;i_phase<2;i_phase++){

      switch(i_phase){
      case 0:
        SID_log("Performing burn-in phase...",SID_LOG_OPEN|SID_LOG_TIMER);
        n_iterations_phase=n_iterations_burn;
        break;
      case 1:
        SID_log("Performing integration phase...",SID_LOG_OPEN|SID_LOG_TIMER);
        n_iterations_phase=n_iterations-n_iterations_burn;
        break;
      }
      
      // Initialize progress reporting
      i_report=0;
      if(n_iterations_phase>20)
        i_iteration_next_report=i_iteration_start+(n_iterations_phase-i_iteration_start)/10;
      else
        i_iteration_next_report=20;

      // Loop until this phase (burn or integrate) is done
      j_covariance =0;
      n_success    =0;
      n_success_all=0;
      n_fail       =0;
      n_fail_all   =0;
      flag_continue=TRUE;
      while(flag_continue){
        // Process one averaging interval at a time
        for(i_avg=0,i_thin=0,i_proposal=0;i_avg<n_avg;){

          // Generate new proposal and determine it's chi^2
          generate_new_MCMC_link(MCMC,P_last,n_P,constant,&RNG,m,b,P_new);
          while(MCMC->map_P_to_M(P_new,MCMC,M_new))
            generate_new_MCMC_link(MCMC,P_last,n_P,constant,&RNG,m,b,P_new);
          MCMC->first_map_call=FALSE;

          i_thin++;
          if(my_chain==SID.My_rank){
            // *** Determine the new proposal's likelihood
            MCMC->compute_MCMC_ln_likelihood(MCMC,M_new,&ln_likelihood_new);
            // *** If we decide to keep this proposal, then ...
            ln_Pr_new=MIN(0.0,ln_likelihood_new-ln_likelihood_last);
            if((double)random_number(&RNG)<=exp(ln_Pr_new)){
              memcpy(P_last,P_new,(size_t)n_P*sizeof(double));      
              for(i_DS=0;i_DS<n_DS;i_DS++)
                memcpy(M_last[i_DS],M_new[i_DS],(size_t)n_M[i_DS]*sizeof(double));
              ln_likelihood_last=ln_likelihood_new;
              ln_Pr_last        =ln_Pr_new;
              flag_success      =TRUE;
            }
            // ... else ...
            else
              flag_success=FALSE;

            // If this proposal is kept by the thining then ...
            if(i_thin>=n_thin){
              i_proposal++;
              if(flag_success)
                n_success++;
              else
                n_fail++;
              if(flag_report_props && my_chain==SID.My_rank){
                if(flag_success)
                  SID_log("Proposal #%09d: SUCCEEDED",SID_LOG_ALLRANKS|SID_LOG_OPEN,i_proposal);
                else
                  SID_log("Proposal #%09d: REJECTED",SID_LOG_ALLRANKS|SID_LOG_OPEN,i_proposal);
                for(i_P=0;i_P<n_P;i_P++)
                  SID_log("%12s=%le (was %le)",SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->P_names[i_P],P_new[i_P],P_last[i_P]);
                SID_log("ln(likelihood) =%le+constant (was %le+constant)",SID_LOG_ALLRANKS|SID_LOG_COMMENT,ln_likelihood_new,ln_likelihood_last);
                SID_log("ln(probability)=%le", SID_LOG_ALLRANKS|SID_LOG_COMMENT,ln_Pr_new);
                SID_log("n_success      =%09d",SID_LOG_ALLRANKS|SID_LOG_COMMENT,n_success);
                SID_log("n_fail         =%09d",SID_LOG_ALLRANKS|SID_LOG_COMMENT,n_fail);
                SID_log("",SID_LOG_ALLRANKS|SID_LOG_NOPRINT|SID_LOG_CLOSE);
              }
              // ... write it to the chain file
              for(i_P=0;i_P<n_P;i_P++){
                P_chain[i_P][i_avg%n_avg]    =P_last[i_P];
                P_proposals[i_P][i_avg%n_avg]=P_new[i_P];
              }
              ln_Pr_chain[i_avg%n_avg]    =ln_Pr_last;
              ln_Pr_proposals[i_avg%n_avg]=ln_Pr_new;
              if(my_chain==SID.My_rank){
                fwrite(&flag_success,sizeof(char),    1,fp_chain);
                fwrite(&ln_Pr_last,  sizeof(double),  1,fp_chain);
                fwrite(P_new,        sizeof(double),n_P,fp_chain);
                for(i_DS=0;i_DS<n_DS;i_DS++)
                  fwrite(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
              }

              // ... and add it to the covariance matrix
              n_covariance++;
              i_covariance++;
              for(i_P=0,k_P=0;i_P<n_P;i_P++){
                P_i_bar_accum[i_P]+=P_last[i_P];
                for(j_P=0;j_P<n_P;j_P++,k_P++)
                  P_ij_bar_accum[k_P]+=P_last[i_P]*P_last[j_P];
              }

              // Update the covariance matrix after each covariance averaging interval
              if(i_covariance>=n_avg_covariance){
                for(i_P=0;i_P<n_P;i_P++)
                  P_i_bar[i_P]=P_i_bar_accum[i_P]/(double)n_covariance;
                for(i_P=0,k_P=0;i_P<n_P;i_P++){
                  for(j_P=0;j_P<n_P;j_P++,k_P++){
                    P_ij_bar[k_P] =P_ij_bar_accum[k_P]/(double)n_covariance;
                    V_compute[k_P]=P_ij_bar[k_P]-P_i_bar[i_P]*P_i_bar[j_P];
                  }
                }
                // Start fresh after the first 2 covariance averaging intervals
                //   of the burn-in phase but keep accumulating otherwise
                if(i_phase==0 && j_covariance<2){
                  for(i_P=0,k_P=0;i_P<n_P;i_P++){
                    P_i_bar_accum[i_P]=0.;
                    P_i_bar[i_P]      =0.;
                    for(j_P=0;j_P<n_P;j_P++,k_P++){
                      P_ij_bar_accum[k_P]=0.;
                      P_ij_bar[k_P]      =0.;
                    }
                  }
                  n_covariance=0;
                }
                // Update MCMC structure during burn-in only
                if(i_phase==0){
                  if(!check_mode_for_flag(MCMC->mode,MCMC_NO_CVM_UPDATE))
                    add_covariance_to_MCMC(MCMC,V_compute);
                }
                i_covariance=0;
                j_covariance++;
              }
            }
          }
          if(i_thin>=n_thin){
            i_thin=0;
            i_avg++;
          }
        } // i_avg
        if(i_phase==0)
          n_iterations_file_burn++;
        n_iterations_file_total++;
        i_iteration++;
        if(flag_report_props && my_chain==SID.My_rank)
          SID_log("Iteration #%07d of #%07d complete for this phase. ",SID_LOG_ALLRANKS|SID_LOG_COMMENT,i_iteration,n_iterations_phase);

        // Generate statistics for the averaging interval we just completed
        if(my_chain==SID.My_rank){
          // Write the statistics to the chain stats file
          compute_MCMC_chain_stats(P_proposals,n_P,n_avg,P_min,P_avg,P_max,dP_avg,auto_cor,slopes,dP_sub,ln_Pr_proposals,&ln_Pr_min,&ln_Pr_avg,&ln_Pr_max);
          fwrite(P_min,   sizeof(double),n_P,fp_stats);
          fwrite(P_avg,   sizeof(double),n_P,fp_stats);
          fwrite(P_max,   sizeof(double),n_P,fp_stats);
          fwrite(dP_avg,  sizeof(double),n_P,fp_stats);
          fwrite(dP_sub,  sizeof(double),n_P,fp_stats);
          fwrite(&ln_Pr_min,sizeof(double),1,fp_stats);
          fwrite(&ln_Pr_avg,sizeof(double),1,fp_stats);
          fwrite(&ln_Pr_max,sizeof(double),1,fp_stats);
          if(flag_autocor_on)
            fwrite(auto_cor,sizeof(double),n_avg-1,fp_stats);
          compute_MCMC_chain_stats(P_chain,n_P,n_avg,P_min,P_avg,P_max,dP_avg,auto_cor,slopes,dP_sub,ln_Pr_chain,&ln_Pr_min,&ln_Pr_avg,&ln_Pr_max);
          fwrite(P_min,   sizeof(double),n_P,fp_stats);
          fwrite(P_avg,   sizeof(double),n_P,fp_stats);
          fwrite(P_max,   sizeof(double),n_P,fp_stats);
          fwrite(dP_avg,  sizeof(double),n_P,fp_stats);
          fwrite(dP_sub,  sizeof(double),n_P,fp_stats);
          fwrite(&ln_Pr_min,sizeof(double),1,fp_stats);
          fwrite(&ln_Pr_avg,sizeof(double),1,fp_stats);
          fwrite(&ln_Pr_max,sizeof(double),1,fp_stats);
          if(flag_autocor_on)
            fwrite(auto_cor,sizeof(double),n_avg-1,fp_stats);
          for(i_P=0;i_P<n_P;i_P++)
            drift[i_P]+=slopes[i_P];
          fwrite(slopes,  sizeof(double),n_P,fp_stats);
          fwrite(drift,   sizeof(double),n_P,fp_stats);
        }
        
        // Check to see if a stop has been called to the run ...
        if(SID.I_am_Master){
          if((fp_stop=fopen(filename_stop,"r"))!=NULL){
            fclose(fp_stop);
            remove(filename_stop);
            flag_stop=TRUE;
          }
        }
        SID_Bcast(&flag_stop,sizeof(int),MASTER_RANK);
        
        // ... if so, stop all ranks and cancel the subsequent analysis stage
        if(flag_stop){
          flag_continue         =FALSE;
          MCMC->flag_analysis_on=FALSE;
          i_phase               =3;
          SID_log("*** Stop file found.  Terminating run. ***",SID_LOG_COMMENT);          
        }
        
        // Check to see if this phase's iterations are complete
        if(i_iteration>=n_iterations_phase)
          flag_continue=FALSE;

        // Report progress
        if(i_iteration==i_iteration_next_report){
          i_report++;
          SID_log("%3d%% complete.",SID_LOG_COMMENT|SID_LOG_TIMER,10*(i_report));
          i_iteration_next_report=MIN(n_iterations_phase,i_iteration_start+(n_iterations_phase-i_iteration_start)*(i_report+1)/10);
        }
      } // while continue
      
      // Report success rate
      SID_Allreduce(&n_success,&n_success_all,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(&n_fail,   &n_fail_all,   1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
      SID_log("Proposal success: %5.3f%% (%lld of %lld)",SID_LOG_COMMENT,1e2*((float)(n_success_all)/(float)(n_success_all+n_fail_all)),n_success_all,n_success_all+n_fail_all);
  
      if(my_chain==SID.My_rank){
        // Write the status of the covariance matrix calculation
        fp_chain_covariance=fopen(filename_chain_covariance,"wb");
        fwrite(&i_covariance, sizeof(int),   1,      fp_chain_covariance);
        fwrite(&n_covariance, sizeof(int),   1,      fp_chain_covariance);
        fwrite(&n_P,          sizeof(int),   1,      fp_chain_covariance);
        fwrite(V_compute,     sizeof(double),n_P*n_P,fp_chain_covariance);
        fwrite(P_i_bar_accum, sizeof(double),n_P,    fp_chain_covariance);
        fwrite(P_ij_bar_accum,sizeof(double),n_P*n_P,fp_chain_covariance);
        fclose(fp_chain_covariance);

        // Reset the covariance matrix calculation if this is the end of the burn-in ...
        if(i_phase==0 && !flag_stop){

          // This is the covariance matrix that will be 
          //   used for the integration phase
          if(!check_mode_for_flag(MCMC->mode,MCMC_NO_CVM_UPDATE))
            add_covariance_to_MCMC(MCMC,V_compute);
          for(i_P=0,k_P=0;i_P<n_P;i_P++){
            P_i_bar_accum[i_P]=0.;
            P_i_bar[i_P]      =0.;
            for(j_P=0;j_P<n_P;j_P++,k_P++){
              P_ij_bar_accum[k_P]=0.;
              P_ij_bar[k_P]      =0.;
            }
          }
          i_covariance=0;
          n_covariance=0;
        }
      }
      i_iteration      =0;
      i_iteration_start=i_iteration;
      SID_log("Done.",SID_LOG_CLOSE);
    }
    if(my_chain==SID.My_rank){
      fp_chain_iterations=fopen(filename_chain_iterations,"wb");
      fwrite(&n_iterations_file_total,sizeof(int),1,fp_chain_iterations);
      fwrite(&n_iterations_file_burn, sizeof(int),1,fp_chain_iterations);
      fclose(fp_chain_iterations);
      fclose(fp_chain);
      fclose(fp_stats);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }

  // Perform analysis on chain(s)
  if(MCMC->flag_analysis_on){
    SID_log("Analyzing chain(s)...",SID_LOG_OPEN);
  
    SID_log("Initializing...",SID_LOG_OPEN);
    for(i_P=0;i_P<n_P;i_P++){
      P_min[i_P]=DBL_MAX;
      P_max[i_P]=DBL_MIN;
      for(j_P=0;j_P<coverage_size;j_P++)
        P_histogram[i_P][j_P]=0;
    }
    for(i_coverage=0;i_coverage<n_coverage;i_coverage++){
      for(j_coverage=0;j_coverage<coverage_size*coverage_size;j_coverage++){
        coverage_true[i_coverage][j_coverage] =0;
        coverage_false[i_coverage][j_coverage]=0;
        coverage_keep[i_coverage][j_coverage] =0;
      }
    }
    for(i_DS=0;i_DS<n_DS;i_DS++){
      for(i_M=0;i_M<n_M[i_DS];i_M++){
        M_min[i_DS][i_M]=DBL_MAX;
        M_max[i_DS][i_M]=DBL_MIN;
        M_avg[i_DS][i_M] =0.;
        dM_avg[i_DS][i_M]=0.;
        for(j_M=0;j_M<coverage_size;j_M++)
          M_histogram[i_DS][i_M][j_M]=0;
      }
    }
    for(i_P=0,k_P=0;i_P<n_P;i_P++){
      P_i_bar_accum[i_P]=0.;
      P_i_bar[i_P]      =0.;
      for(j_P=0;j_P<n_P;j_P++,k_P++){
        P_ij_bar_accum[k_P]=0.;
        P_ij_bar[k_P]      =0.;
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Process the chain(s)
    SID_log("Process chain file(s)...",SID_LOG_OPEN);
    if(my_chain==SID.My_rank){
      fp_chain=fopen(filename_chain,"rb");
      for(i_iteration=0,flag_init=TRUE,n_used=0;i_iteration<n_iterations;i_iteration++){
        for(i_avg=0;i_avg<n_avg;i_avg++){
          fread(&flag_success,sizeof(char),  1,  fp_chain);
          fread(&ln_Pr_last,  sizeof(double),1,  fp_chain);
          fread(P_new,        sizeof(double),n_P,fp_chain);
          for(i_DS=0;i_DS<n_DS;i_DS++)
            fread(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
          if(i_iteration>=n_iterations_burn){
            if(flag_init || flag_success){
              memcpy(P_last,P_new,(size_t)n_P*sizeof(double));      
              for(i_DS=0;i_DS<n_DS;i_DS++)
                memcpy(M_last[i_DS],M_new[i_DS],(size_t)n_M[i_DS]*sizeof(double));
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
      // Finish averages
      for(i_P=0;i_P<n_P;i_P++)
        P_avg[i_P]/=(double)n_used;
      for(i_DS=0;i_DS<n_DS;i_DS++){
        for(i_P=0;i_P<n_M[i_DS];i_P++)
          M_avg[i_DS][i_P]/=(double)n_used;
      }
      rewind(fp_chain);
      for(i_iteration=0,flag_init=TRUE;i_iteration<n_iterations;i_iteration++){
        for(i_avg=0;i_avg<n_avg;i_avg++){
          fread(&flag_success,sizeof(char),  1,  fp_chain);
          fread(&ln_Pr_last,  sizeof(double),1,  fp_chain);
          fread(P_new,        sizeof(double),n_P,fp_chain);
          for(i_DS=0;i_DS<n_DS;i_DS++)
            fread(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
          if(i_iteration>=n_iterations_burn){
            if(flag_init || flag_success){
              memcpy(P_last,P_new,(size_t)n_P*sizeof(double));      
              for(i_DS=0;i_DS<n_DS;i_DS++)
                memcpy(M_last[i_DS],M_new[i_DS],(size_t)n_M[i_DS]*sizeof(double));
              flag_init=FALSE;
            }

            // Build coverage maps
            for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
              bin_x=(double)(coverage_size)*(P_new[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
              for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
                bin_y=(double)(coverage_size)*(P_new[j_P]-P_min[j_P])/(P_max[j_P]-P_min[j_P]);
                if(bin_x>=0 && bin_x<coverage_size && bin_y>=0 && bin_y<coverage_size){
                  switch(flag_success){
                  case TRUE:
                    coverage_true[i_coverage][bin_x*coverage_size+bin_y]++;
                    break;
                  case FALSE:
                    coverage_false[i_coverage][bin_x*coverage_size+bin_y]++;
                    break;
                  default:
                    SID_trap_error("Unknown success flag (%d) when constructing coverage map.",ERROR_LOGIC,(int)flag_success);
                    break;
                  }
                }
              }
            }
            for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
              bin_x=(double)(coverage_size)*(P_last[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
              for(j_P=i_P+1;j_P<n_P;j_P++,i_coverage++){
                bin_y=(double)(coverage_size)*(P_last[j_P]-P_min[j_P])/(P_max[j_P]-P_min[j_P]);
                if(bin_x>=0 && bin_x<coverage_size && bin_y>=0 && bin_y<coverage_size)
                  coverage_keep[i_coverage][bin_x*coverage_size+bin_y]++;
              }
            }

            // Build histograms
            for(i_P=0,i_coverage=0;i_P<n_P;i_P++){
              bin_x=(double)(coverage_size)*(P_last[i_P]-P_min[i_P])/(P_max[i_P]-P_min[i_P]);
              if(bin_x>=0 && bin_x<coverage_size)
                P_histogram[i_P][bin_x]++;
            }
            for(i_DS=0;i_DS<n_DS;i_DS++){
              for(i_M=0;i_M<n_M[i_DS];i_M++){
                bin_x=(double)(coverage_size)*(M_last[i_DS][i_M]-M_min[i_DS][i_M])/(M_max[i_DS][i_M]-M_min[i_DS][i_M]);
                if(bin_x>=0 && bin_x<coverage_size)
                  M_histogram[i_DS][i_M][bin_x]++;
              }
            }
      
            // Compute parameter standard deviations
            for(i_P=0;i_P<n_P;i_P++)
              dP_avg[i_P]+=pow(P_last[i_P]-P_avg[i_P],2.);
            for(i_DS=0;i_DS<n_DS;i_DS++){
              for(i_P=0;i_P<n_M[i_DS];i_P++)
                dM_avg[i_DS][i_P]+=pow(M_last[i_DS][i_P]-M_avg[i_DS][i_P],2.);
            }
          }
        }
      } // i_iteration
      fclose(fp_chain);
      // Finish standard deviations
      for(i_DS=0;i_DS<n_DS;i_DS++){
        for(i_P=0;i_P<n_P;i_P++){
          if(i_DS==0)
            dP_avg[i_P]=sqrt(dP_avg[i_P]/(double)n_used);
          dM_avg[i_DS][i_P]=sqrt(dM_avg[i_DS][i_P]/(double)n_used);
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

      // Compute parameter confidence intervals from histograms
      for(i_P=0;i_P<n_P;i_P++){
        merge_sort(P_histogram[i_P],(size_t)coverage_size,&histogram_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
        accumulator       =P_histogram[i_P][histogram_index[coverage_size-1]];
        P_best_index      =histogram_index[coverage_size-1];
        P_lo_68_index     =histogram_index[coverage_size-1];
        P_hi_68_index     =histogram_index[coverage_size-1];
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
          fwrite(&P_contour_68[i_coverage], sizeof(double),1,                          fp_coverage);
          fwrite(&P_contour_95[i_coverage], sizeof(double),1,                          fp_coverage);
          fwrite(coverage_true[i_coverage], sizeof(size_t),coverage_size*coverage_size,fp_coverage);
          fwrite(coverage_false[i_coverage],sizeof(size_t),coverage_size*coverage_size,fp_coverage);
          fwrite(coverage_keep[i_coverage], sizeof(size_t),coverage_size*coverage_size,fp_coverage);
        }
      }
      fclose(fp_coverage);
      SID_log("Done.",SID_LOG_CLOSE);

      // Write coverage maps
      SID_log("Writing histograms...",SID_LOG_OPEN);
      fp_histograms=fopen(filename_histograms,"wb");
      for(i_P=0;i_P<n_P;i_P++){
        fwrite(&P_best[i_P],    sizeof(double),1,            fp_histograms);
        fwrite(&P_lo_68[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_hi_68[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_lo_95[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(&P_hi_95[i_P],   sizeof(double),1,            fp_histograms);
        fwrite(P_histogram[i_P],sizeof(size_t),coverage_size,fp_histograms);
      }
      for(i_DS=0;i_DS<n_DS;i_DS++){
        for(i_M=0;i_M<n_M[i_DS];i_M++){
          fwrite(&M_best[i_DS][i_M],    sizeof(double),1,            fp_histograms);
          fwrite(&M_lo_68[i_DS][i_M],   sizeof(double),1,            fp_histograms);
          fwrite(&M_hi_68[i_DS][i_M],   sizeof(double),1,            fp_histograms);
          fwrite(&M_lo_95[i_DS][i_M],   sizeof(double),1,            fp_histograms);
          fwrite(&M_hi_95[i_DS][i_M],   sizeof(double),1,            fp_histograms);
          fwrite(M_histogram[i_DS][i_M],sizeof(size_t),coverage_size,fp_histograms);
        }
      }
      fclose(fp_histograms);
      SID_log("Done.",SID_LOG_CLOSE);

      // Write results
      SID_log("Writing final statistics...",SID_LOG_OPEN);
      sprintf(filename_results,"%s/fit_for_parameters.dat",filename_results_dir);
      fp_results=fopen(filename_results,"w");
      fprintf(fp_results,"# MCMC parameter fit results to %s\n",MCMC->problem_name);
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
      for(i_P=0;i_P<n_P;i_P++){
        fprintf(fp_results,"%s",MCMC->P_names[i_P]);
        for(i_array=0;i_array<MCMC->n_arrays;i_array++)
          fprintf(fp_results," %le",MCMC->array[i_array][i_P]);
        fprintf(fp_results," %le %le %le %le %le %le %le %le\n",MCMC->P_init[i_P],P_avg[i_P],dP_avg[i_P],P_best[i_P],P_lo_68[i_P],P_hi_68[i_P],P_lo_95[i_P],P_hi_95[i_P]);
      }
      fclose(fp_results);
      current_DS=MCMC->DS;
      i_DS      =0;
      while(current_DS!=NULL){
        next_DS  =current_DS->next;
        M_target =current_DS->M_target;
        dM_target=current_DS->dM_target;
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
        fprintf(fp_results,"# %s (%02d) Average\n",                      column_txt,i_column++);
        fprintf(fp_results,"# %s (%02d) Standard deviation\n",           column_txt,i_column++);
        fprintf(fp_results,"# %s (%02d) Maximum likelihood value\n",     column_txt,i_column++);
        fprintf(fp_results,"# %s (%02d) Lower limit (68%% confidence)\n",column_txt,i_column++);
        fprintf(fp_results,"# %s (%02d) Upper limit (68%% confidence)\n",column_txt,i_column++);
        fprintf(fp_results,"# %s (%02d) Lower limit (95%% confidence)\n",column_txt,i_column++);
        fprintf(fp_results,"# %s (%02d) Upper limit (95%% confidence)\n",column_txt,i_column++);
        for(i_P=0;i_P<n_M[i_DS];i_P++){
          for(i_array=0;i_array<n_M_arrays[i_DS];i_array++)
            fprintf(fp_results,"%le ",M_arrays[i_DS][i_array][i_P]);
          fprintf(fp_results,"%le %le %le %le %le %le %le %le %le\n",M_target[i_P],dM_target[i_P],M_avg[i_DS][i_P],dM_avg[i_DS][i_P],M_best[i_DS][i_P],M_lo_68[i_DS][i_P],M_hi_68[i_DS][i_P],M_lo_95[i_DS][i_P],M_hi_95[i_DS][i_P]);
        }
        fclose(fp_results);    
        current_DS=next_DS;
        i_DS++;
      }    
      SID_log("Done.",SID_LOG_CLOSE);
      SID_log("Done.",SID_LOG_CLOSE);
    }
  }

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  free_RNG(&RNG);
  gsl_vector_free(b);
  SID_free(SID_FARG P_i_bar_accum);
  SID_free(SID_FARG P_ij_bar_accum);
  SID_free(SID_FARG P_i_bar);
  SID_free(SID_FARG P_ij_bar);
  SID_free(SID_FARG V_compute);
  for(i_P=0;i_P<n_coverage;i_P++){
    SID_free(SID_FARG coverage_true[i_P]);
    SID_free(SID_FARG coverage_false[i_P]);
    SID_free(SID_FARG coverage_keep[i_P]);
  }
  for(i_P=0;i_P<n_P;i_P++)
    SID_free(SID_FARG P_histogram[i_P]);
  SID_free(SID_FARG P_histogram);
  SID_free(SID_FARG P_contour_68);
  SID_free(SID_FARG P_contour_95);
  SID_free(SID_FARG P_best);
  SID_free(SID_FARG coverage_true);
  SID_free(SID_FARG coverage_false);
  SID_free(SID_FARG coverage_keep);
  for(i_DS=0;i_DS<n_DS;i_DS++){
    SID_free(SID_FARG M_new[i_DS]);
    SID_free(SID_FARG M_last[i_DS]);
    SID_free(SID_FARG M_avg[i_DS]);
    SID_free(SID_FARG dM_avg[i_DS]);
    SID_free(SID_FARG M_best[i_DS]);
    SID_free(SID_FARG M_min[i_DS]);
    SID_free(SID_FARG M_max[i_DS]);
    SID_free(SID_FARG M_lo_68[i_DS]);
    SID_free(SID_FARG M_hi_68[i_DS]);
    SID_free(SID_FARG M_lo_95[i_DS]);
    SID_free(SID_FARG M_hi_95[i_DS]);
  }
  SID_free(SID_FARG n_M);
  SID_free(SID_FARG M_new);
  SID_free(SID_FARG M_last);
  SID_free(SID_FARG M_avg);
  SID_free(SID_FARG dM_avg);
  SID_free(SID_FARG M_best);
  SID_free(SID_FARG M_min);
  SID_free(SID_FARG M_max);
  SID_free(SID_FARG M_lo_68);
  SID_free(SID_FARG M_hi_68);
  SID_free(SID_FARG M_lo_95);
  SID_free(SID_FARG M_hi_95);
  SID_free(SID_FARG P_lo_68);
  SID_free(SID_FARG P_hi_68);
  SID_free(SID_FARG P_lo_95);
  SID_free(SID_FARG P_hi_95);
  SID_free(SID_FARG P_min);
  SID_free(SID_FARG P_max);
  SID_free(SID_FARG P_new);
  SID_free(SID_FARG P_last);
  SID_free(SID_FARG P_avg);
  SID_free(SID_FARG dP_avg);
  SID_free(SID_FARG slopes);
  SID_free(SID_FARG dP_sub);
  SID_free(SID_FARG drift);
  for(i_P=0;i_P<n_P;i_P++){
    SID_free(SID_FARG P_chain[i_P]);
    SID_free(SID_FARG P_proposals[i_P]);
    if(flag_autocor_on)
      SID_free(SID_FARG auto_cor[i_P]);
  }
  SID_free(SID_FARG P_chain);
  SID_free(SID_FARG ln_Pr_chain);
  SID_free(SID_FARG ln_Pr_proposals);
  SID_free(SID_FARG P_proposals);
  SID_free(SID_FARG n_M_arrays);
  SID_free(SID_FARG M_arrays);
  if(flag_autocor_on)
    SID_free(SID_FARG auto_cor);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}
