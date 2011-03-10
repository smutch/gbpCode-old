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

//todo: get rid of constant
//      write T to a file somewhere for restarts
//      deal with covariance in restarts

void compute_MCMC(MCMC_info *MCMC){
  char      filename_output_dir[MAX_FILENAME_LENGTH];
  char      filename_chain_dir[MAX_FILENAME_LENGTH];
  char      filename_results_dir[MAX_FILENAME_LENGTH];
  char      filename_plots_dir[MAX_FILENAME_LENGTH];
  char      filename_run[MAX_FILENAME_LENGTH];
  char      filename_chain[MAX_FILENAME_LENGTH];
  char      filename_chain_config[MAX_FILENAME_LENGTH];
  char      filename_chain_covariance[MAX_FILENAME_LENGTH];
  char      filename_stats[MAX_FILENAME_LENGTH];
  char      filename_coverage[MAX_FILENAME_LENGTH];
  char      filename_histograms[MAX_FILENAME_LENGTH];
  char      filename_results[MAX_FILENAME_LENGTH];
  char      filename_stop[MAX_FILENAME_LENGTH];
  char      column_txt[MAX_FILENAME_LENGTH];
  char      problem_name_test[MCMC_NAME_SIZE];
  char      P_name_test[MCMC_NAME_SIZE];
  char      format_string[64];
  char      name_test[MCMC_NAME_SIZE];
  char      array_name_test[MCMC_NAME_SIZE];
  int       n_avg_test;
  int       flag_autocor_on_test;
  int       n_P_test;
  int       n_DS_test;
  int       n_arrays_test;
  int       n_M_test;
  double    P_init_test;
  double    M_target_test;
  double    dM_target_test;
  double    array_test;
  double    P_min_test;
  double    P_max_test;
  double    ln_Pr_min;
  double    ln_Pr_max;
  double    ln_Pr_avg;
  int       i_tune;
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
  int       bin_x;
  int       bin_y;
  int       n_avg;
  int       n_avg_burn;
  int       n_avg_integrate;
  int       n_P,n_C;
  int       n_gals;
  double   *x_P;
  double  **x_M;
  double   *M_target;
  double   *dM_target;
  double   *P_min;
  double   *P_max;
  double   *P_new;
  double   *P_last;
  double   *P_init;
  double   *P_chain;
  double   *P_avg;
  double   *dP_avg;
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
  double  **P_chain_stats;
  double  **P_prop_stats;
  size_t **P_histogram;
  size_t **coverage_true;
  size_t **coverage_false;
  size_t **coverage_keep;
  double  *V_read;
  double    L_x,L_y,L_z;
  int       n_x,n_y,n_z;
  int       n_C_x;
  FILE     *fp_run;
  FILE     *fp_chain;
  FILE     *fp_chain_config;
  FILE     *fp_chain_covariance;
  FILE     *fp_stats;
  FILE     *fp_coverage;
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
  int       flag_autocor_on;
  int       n_coverage;
  int       coverage_size;
  int       i_phase;
  int       i_array;
  int      *n_M;
  int       flag_initialized;
  int       flag_no_map_write;
  int       flag_no_map_write_test;
  int       n_used;
  double      RN;
  double       ***P_arrays;
  int            *n_P_arrays;
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
  int             n_iterations_file_total;
  int             n_iterations_file_burn;
  int             flag_restart=FALSE;
  MCMC_DS_info   *current_DS;
  MCMC_DS_info   *next_DS;

  SID_log("Performing MCMC...",SID_LOG_OPEN|SID_LOG_TIMER);

  SID_log("Initializing...",SID_LOG_OPEN);

  // Initialize dataset arrays
  init_MCMC_DS(MCMC);

  // Grab some stuff from the MCMC structure
  n_P                   =MCMC->n_P;
  n_DS                  =MCMC->n_DS;
  n_avg                 =MCMC->n_avg;
  n_thin                =MCMC->n_thin;
  coverage_size         =MCMC->coverage_size;
  flag_autocor_on       =MCMC->flag_autocor_on;
  my_chain              =MCMC->my_chain;
  n_iterations          =MCMC->n_iterations;
  n_iterations_burn     =MCMC->n_iterations_burn;
  n_iterations_integrate=n_iterations-n_iterations_burn;
  flag_no_map_write     =check_mode_for_flag(MCMC->mode,MCMC_MODE_NO_MAP_WRITE);

  SID_log("temperature            = %lf",SID_LOG_COMMENT,MCMC->temperature);
  SID_log("n_avg                  = %d",SID_LOG_COMMENT,n_avg);
  SID_log("n_thin                 = %d",SID_LOG_COMMENT,n_thin);
  SID_log("coverage_size          = %d",SID_LOG_COMMENT,coverage_size);
  SID_log("n_iterations           = %d",SID_LOG_COMMENT,n_iterations);
  SID_log("n_iterations_burn      = %d",SID_LOG_COMMENT,n_iterations_burn);
  SID_log("n_iterations_integrate = %d",SID_LOG_COMMENT,n_iterations_integrate);

  // Initialize random number generator
  SID_log("random seed            = %d",SID_LOG_COMMENT,MCMC->seed);

  if(MCMC->V!=NULL)
    SID_log("Covariance matrix is initialized.",SID_LOG_COMMENT);
  else
    SID_log("Covariance matrix is NOT initialized.",SID_LOG_COMMENT);
  if(flag_autocor_on)
    SID_log("Auto-correlation  is on.",SID_LOG_COMMENT);
  else  
    SID_log("Auto-correlation  is off.",SID_LOG_COMMENT);

  // Initialize chain arrays
  n_M            =MCMC->n_M;
  M_new          =MCMC->M_new;
  M_last         =MCMC->M_last;
  P_new          =MCMC->P_new;
  P_last         =MCMC->P_last;
  P_init         =MCMC->P_init;
  P_chain        =MCMC->P_chain;
  P_min          =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_max          =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_avg          =(double  *)SID_malloc(sizeof(double)  *n_P);
  dP_avg         =(double  *)SID_malloc(sizeof(double)  *n_P);
  ln_Pr_chain    =(double  *)SID_malloc(sizeof(double)  *n_avg);
  ln_Pr_proposals=(double  *)SID_malloc(sizeof(double)  *n_avg);
  drift          =(double  *)SID_malloc(sizeof(double)  *n_P);
  slopes         =(double  *)SID_malloc(sizeof(double)  *n_P);
  dP_sub         =(double  *)SID_malloc(sizeof(double)  *n_P);
  P_chain_stats  =(double **)SID_malloc(sizeof(double *)*n_P);
  P_prop_stats   =(double **)SID_malloc(sizeof(double *)*n_P);
  if(flag_autocor_on)
    auto_cor=(double **)SID_malloc(sizeof(double *)*n_P);
  else
    auto_cor=NULL;
  for(i_P=0;i_P<n_P;i_P++){
    drift[i_P]        =0.;
    slopes[i_P]       =0.;
    dP_sub[i_P]       =0.;
    P_chain_stats[i_P]=(double *)SID_malloc(sizeof(double)*n_avg);
    P_prop_stats[i_P] =(double *)SID_malloc(sizeof(double)*n_avg);
    if(flag_autocor_on)
      auto_cor[i_P]=(double *)SID_malloc(sizeof(double)*n_avg);
  }

  // Set directories
  sprintf(filename_output_dir, "%s/",        MCMC->filename_output_dir);
  sprintf(filename_chain_dir,  "%s/chains/", filename_output_dir);
  sprintf(filename_results_dir,"%s/results/",filename_output_dir);
  sprintf(filename_plots_dir,  "%s/plots/",  filename_output_dir);

  // Set filenames
  sprintf(filename_run,             "%s/run.dat",                  filename_output_dir);
  sprintf(filename_chain,           "%s/chain_trace_%06d.dat",     filename_chain_dir,my_chain);
  sprintf(filename_chain_config,    "%s/chain_config_%06d.dat",    filename_chain_dir,my_chain);
  sprintf(filename_chain_covariance,"%s/chain_covariance_%06d.dat",filename_chain_dir,my_chain);
  sprintf(filename_stats,           "%s/chain_stats_%06d.dat",     filename_chain_dir,my_chain);
  sprintf(filename_coverage,        "%s/coverage.dat",             filename_results_dir);
  sprintf(filename_histograms,      "%s/histograms.dat",           filename_results_dir);
  sprintf(filename_stop,            "%s/stop",                     filename_output_dir);
  
  SID_log("Done.",SID_LOG_CLOSE);

  // Perform integration
    SID_log("Initializing integration...",SID_LOG_OPEN);
    MCMC->flag_integrate_on=TRUE;    

    // Make sure the needed directories exist
    mkdir(filename_output_dir, 02755);
    mkdir(filename_chain_dir,  02755);
    mkdir(filename_results_dir,02755);
    mkdir(filename_plots_dir,  02755);

    // Read/Write Header file
    if(SID.I_am_Master){
      // If run.dat already exists, this is a restart.  Check that it is consistant with the current run.
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
        fread(&flag_autocor_on_test,sizeof(int),1,fp_run);
        if(flag_autocor_on_test!=flag_autocor_on)
          SID_trap_error("Autocorrelation flags are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,flag_autocor_on,flag_autocor_on_test);
        fread(&flag_no_map_write_test,sizeof(int),1,fp_run);
        if(flag_no_map_write_test!=flag_no_map_write)
          SID_trap_error("Map write flags are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,flag_no_map_write,flag_no_map_write_test);
        fread(&n_P_test,sizeof(int),1,fp_run);
        if(n_P_test!=n_P)
          SID_trap_error("The number of paramaters is inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_P,n_P_test);
        MCMC->P_name_length=0;
        for(i_P=0;i_P<n_P;i_P++){
          fread(P_name_test, sizeof(char),MCMC_NAME_SIZE,fp_run);
          if(strcmp(P_name_test,MCMC->P_names[i_P]))
            SID_trap_error("Parameter #%d's names are inconsistant (i.e. {%s}!={%s}).",ERROR_LOGIC,i_P,MCMC->P_names[i_P],P_name_test);
          fread(&P_init_test,sizeof(double),1,fp_run);
          if(P_init_test!=MCMC->P_init[i_P])
            SID_trap_error("Parameter #%d's initial values are inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_P,MCMC->P_init[i_P],P_init_test);
          fread(&P_min_test,sizeof(double),1,fp_run);
          if(P_min_test!=MCMC->P_limit_min[i_P])
            SID_trap_error("Parameter #%d's minimum values are inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_P,MCMC->P_limit_min[i_P],P_min_test);
          fread(&P_max_test,sizeof(double),1,fp_run);
          if(P_max_test!=MCMC->P_limit_max[i_P])
            SID_trap_error("Parameter #%d's maximum values are inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_P,MCMC->P_limit_max[i_P],P_max_test);
          MCMC->P_name_length=MAX(MCMC->P_name_length,strlen(MCMC->P_names[i_P]));
        }
        sprintf(MCMC->P_name_format,"%%-%ds",MCMC->P_name_length);
        fread(&n_arrays_test,sizeof(int),1,fp_run);
        if(n_arrays_test!=MCMC->n_arrays)
          SID_trap_error("Numbers of project arrays are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,MCMC->n_arrays,n_arrays_test);
        for(i_array=0;i_array<MCMC->n_arrays;i_array++){
          fread(&array_name_test,sizeof(char),MCMC_NAME_SIZE,fp_run);
          if(strcmp(array_name_test,MCMC->array_name[i_array]))
            SID_trap_error("Project array names are inconsisitant (i.e. {%s}!={%s}).",ERROR_LOGIC,MCMC->array_name[i_array],array_name_test);
          for(i_P=0;i_P<n_P;i_P++){
            fread(&array_test,sizeof(double),1,fp_run);
            if(array_test!=MCMC->array[i_array][i_P])
              SID_trap_error("Project array #%d element #%d is inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_array,i_P,MCMC->array[i_array][i_P],array_test);
          }
        }
        fread(&n_DS_test,sizeof(int),1,fp_run);
        if(n_DS_test!=n_DS)
          SID_trap_error("The number of datasets is inconsistant (i.e. %d!=%d).",ERROR_LOGIC,n_DS,n_DS_test);
        current_DS=MCMC->DS;
        i_DS=0;
        while(current_DS!=NULL){
          next_DS=current_DS->next;
          fread(name_test,     sizeof(char),  MCMC_NAME_SIZE, fp_run);
          if(strcmp(name_test,current_DS->name))
            SID_trap_error("Dataset #%d's names are inconsistant (i.e. {%s}!={%s}).",ERROR_LOGIC,i_DS,current_DS->name,name_test);
          fread(&n_M_test,sizeof(int),1,fp_run);
          if(n_M_test!=n_M[i_DS])
            SID_trap_error("The sizes of dataset #%d are inconsistant (i.e. %d!=%d).",ERROR_LOGIC,i_DS,n_M_test,n_M[i_DS]);
          for(i_M=0;i_M<current_DS->n_M;i_M++){
            fread(&M_target_test, sizeof(double),1,fp_run);
            if(M_target_test!=current_DS->M_target[i_M])
              SID_trap_error("Dataset #%d, element #%d is inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_DS,i_M,current_DS->M_target[i_M],M_target_test);
          }
          for(i_M=0;i_M<current_DS->n_M;i_M++){
            fread(&dM_target_test,sizeof(double),1,fp_run);
            if(dM_target_test!=current_DS->dM_target[i_M])
              SID_trap_error("Dataset #%d, uncertainty element #%d is inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_DS,i_M,current_DS->M_target[i_M],M_target_test);
          }
          fread(&n_arrays_test,sizeof(int),1,fp_run);
          if(n_arrays_test!=current_DS->n_arrays)
            SID_trap_error("The number of arrays in dataset #%d is inconsistant (i.e. %d!=%d).",ERROR_LOGIC,i_DS,n_arrays_test,current_DS->n_arrays);
          for(i_array=0;i_array<current_DS->n_arrays;i_array++){
            fread(array_name_test,sizeof(char),MCMC_NAME_SIZE,fp_run);
            if(strcmp(array_name_test,current_DS->array_name[i_array]))
              SID_trap_error("Array name #%d for dataset #%d is inconsisitant (i.e. {%s}!={%s}).",ERROR_LOGIC,i_array,i_DS,current_DS->array_name[i_array],array_name_test);
            for(i_M=0;i_M<current_DS->n_M;i_M++){
              fread(&array_test,sizeof(double),1,fp_run);
              if(array_test!=current_DS->array[i_array][i_M])
                SID_trap_error("Array #%d, element #%d for dataset #%d is inconsistant (i.e. %le!=%le).",ERROR_LOGIC,i_array,i_M,i_DS,current_DS->array[i_array][i_M],array_test);
            }
          }
          current_DS=next_DS;
          i_DS++;
        }
        fclose(fp_run);
        SID_log("Done.",SID_LOG_CLOSE);
      }
      // ... else if run.dat does not exist, create it
      else{
        SID_log("Write header file...",SID_LOG_OPEN);
        fp_run=fopen(filename_run,"wb");
        // Stuff relating to this MCMC project
        fwrite(MCMC->problem_name,sizeof(char),MCMC_NAME_SIZE,fp_run);
        fwrite(&n_avg,            sizeof(int),   1,   fp_run);
        fwrite(&flag_autocor_on,  sizeof(int),   1,   fp_run);
        fwrite(&flag_no_map_write,sizeof(int),   1,   fp_run);
        fwrite(&n_P,              sizeof(int),   1,   fp_run);
        for(i_P=0;i_P<n_P;i_P++){
          fwrite(MCMC->P_names[i_P],  sizeof(char),  MCMC_NAME_SIZE,fp_run);
          fwrite(&(MCMC->P_init[i_P]),sizeof(double),1,             fp_run);
          fwrite(&(MCMC->P_limit_min[i_P]), sizeof(double),1,             fp_run);
          fwrite(&(MCMC->P_limit_max[i_P]), sizeof(double),1,             fp_run);
        }
        fwrite(&(MCMC->n_arrays),sizeof(int),1,fp_run);
        for(i_array=0;i_array<MCMC->n_arrays;i_array++){
          fwrite(MCMC->array_name[i_array],sizeof(char),  MCMC_NAME_SIZE,fp_run);
          fwrite(MCMC->array[i_array],     sizeof(double),n_P,           fp_run);
        }
        // Stuff relating to the constraining datasets
        fwrite(&n_DS,sizeof(int),1,   fp_run);
        current_DS=MCMC->DS;
        while(current_DS!=NULL){
          next_DS=current_DS->next;
          fwrite(current_DS->name,       sizeof(char),   MCMC_NAME_SIZE,fp_run);
          fwrite(&(current_DS->n_M),     sizeof(int),                 1,fp_run);
          fwrite(current_DS->M_target,   sizeof(double),current_DS->n_M,fp_run);
          fwrite(current_DS->dM_target,  sizeof(double),current_DS->n_M,fp_run);
          fwrite(&(current_DS->n_arrays),sizeof(int),                 1,fp_run);
          for(i_array=0;i_array<current_DS->n_arrays;i_array++){
            fwrite(current_DS->array_name[i_array],sizeof(char),  MCMC_NAME_SIZE, fp_run);
            fwrite(current_DS->array[i_array],     sizeof(double),current_DS->n_M,fp_run);
          }
          current_DS=next_DS;
        }
        fclose(fp_run);
        SID_log("Done.",SID_LOG_CLOSE);
      }
    }
    SID_Bcast(&flag_restart,sizeof(int),MASTER_RANK,SID.COMM_WORLD);

    // If this is NOT a restart, start from scratch ...
    n_iterations_file_total=0;
    n_iterations_file_burn =0;
    if(!flag_restart){
      SID_log("Initializing state...",SID_LOG_OPEN);

      // Report the starting conditions 
      if(my_chain==SID.My_rank){
        SID_log("Initial parameters:",SID_LOG_ALLRANKS|SID_LOG_OPEN);
        sprintf(format_string,"%s = %%13.6le",MCMC->P_name_format);
        for(i_P=0;i_P<n_P;i_P++)
          SID_log(format_string,SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->P_names[i_P],P_last[i_P]);
        SID_log("",SID_LOG_ALLRANKS|SID_LOG_NOPRINT|SID_LOG_CLOSE);
      }

      // Perform autotuning (if requested)
      if(check_mode_for_flag(MCMC->mode,MCMC_MODE_AUTOTUNE))
        autotune_MCMC(MCMC);
      MCMC->flag_init_chain=TRUE;

      // Set the initial state
      SID_log("Writing chain config file...",SID_LOG_OPEN);
      fp_chain             =fopen(filename_chain,"wb");
      fp_stats             =fopen(filename_stats,"wb");
      fp_chain_config      =fopen(filename_chain_config,"wb");
      fwrite(&n_iterations_file_total,sizeof(int),   1,      fp_chain_config);
      fwrite(&n_iterations_file_burn, sizeof(int),   1,      fp_chain_config);
      fwrite(&(MCMC->temperature),    sizeof(double),1,      fp_chain_config);
      fwrite(MCMC->V,                 sizeof(double),n_P*n_P,fp_chain_config);
      fclose(fp_chain_config);
      SID_log("Done.",SID_LOG_CLOSE);

      SID_log("Done.",SID_LOG_CLOSE);
    }
    // ... else read the state we left-off from ...
    else if(my_chain==SID.My_rank){
      SID_log("Loading previous state...",SID_LOG_OPEN);
      MCMC->flag_init_chain=FALSE;
      // ... fetch the number of intervals that have already been computed ...
      n_iterations_file_total=0;
      if((fp_chain_config=fopen(filename_chain_config,"rb"))!=NULL){
        SID_log("Reading the existant number of iterations...",SID_LOG_OPEN);
        V_read=(double *)SID_malloc(sizeof(double)*n_P*n_P);
        fread(&n_iterations_file_total,sizeof(int),   1,      fp_chain_config);
        fread(&n_iterations_file_burn, sizeof(int),   1,      fp_chain_config);
        fread(&(MCMC->temperature),    sizeof(double),1,      fp_chain_config);
        fread(V_read,                  sizeof(double),n_P*n_P,fp_chain_config);
        set_MCMC_covariance(MCMC,V_read);
        SID_free(SID_FARG V_read);
        SID_log("# burn  iterations = %d (%d requested)",SID_LOG_COMMENT,n_iterations_file_burn ,n_iterations_burn);
        SID_log("# total iterations = %d (%d requested)",SID_LOG_COMMENT,n_iterations_file_total,n_iterations);
        SID_log("Temperature        = %le",              SID_LOG_COMMENT,MCMC->temperature);
        fclose(fp_chain_config);
        SID_log("Done.",SID_LOG_CLOSE);
      }

      // ... scan to the end of the chain and read the last-used parameter set
      if(n_iterations_file_total>0){
        SID_log("Reading the last-used parameter set...",SID_LOG_OPEN);
        fp_chain=fopen(filename_chain,"rb");
        for(i_iteration=0;i_iteration<n_iterations_file_total;i_iteration++){
          for(i_avg=0;i_avg<n_avg;i_avg++){
            fread(&flag_success,sizeof(char),  1,fp_chain);
            if((i_iteration+i_avg)!=0){
              MCMC->ln_likelihood_last=MCMC->ln_likelihood_new;
              memcpy(P_last,P_new,n_P*sizeof(double));
              if(!flag_no_map_write){
                for(i_DS=0;i_DS<n_DS;i_DS++)
                  memcpy(M_last[i_DS],M_new[i_DS],n_M[i_DS]*sizeof(double));
              }
            }
            fread(&(MCMC->ln_likelihood_new),sizeof(double),1,fp_chain);
            fread(P_new,sizeof(double),n_P,fp_chain);
            if(!flag_no_map_write){
              for(i_DS=0;i_DS<n_DS;i_DS++)
                fread(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
            }
            if((i_iteration+i_avg)==0){
              MCMC->ln_likelihood_chain=MCMC->ln_likelihood_new;
              MCMC->ln_likelihood_last =MCMC->ln_likelihood_new;
            }
            if(flag_success){
              MCMC->ln_likelihood_chain=MCMC->ln_likelihood_new;
              memcpy(P_chain,P_new,n_P*sizeof(double));
            }
          }        
        }
        fclose(fp_chain);
        fp_chain=fopen(filename_chain,"ab");
        SID_log("Done.",SID_LOG_CLOSE);

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
        fp_stats=fopen(filename_stats,"ab");        
        SID_log("Done.",SID_LOG_CLOSE);
      }
      else{
        memcpy(P_last,P_init,n_P*sizeof(double));
        MCMC->flag_init_chain=TRUE;
        fp_chain=fopen(filename_chain,"wb");
        fp_stats=fopen(filename_stats,"wb");        
      }

      // Report the starting conditions 
      if(my_chain==SID.My_rank){
        SID_log("Resume from parameters:",SID_LOG_ALLRANKS|SID_LOG_OPEN);
        sprintf(format_string,"%s = %%13.6le",MCMC->P_name_format);
        for(i_P=0;i_P<n_P;i_P++)
          SID_log(format_string,SID_LOG_ALLRANKS|SID_LOG_COMMENT,MCMC->P_names[i_P],P_last[i_P]);
        SID_log("",SID_LOG_ALLRANKS|SID_LOG_NOPRINT|SID_LOG_CLOSE);
      }

      SID_log("Done.",SID_LOG_CLOSE);

    } // End of restart stuff

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

    // This is the end of the initialization
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
      flag_continue=TRUE;
      while(flag_continue){
        // Process one averaging interval at a time
        for(i_avg=0,i_thin=1;i_avg<n_avg;i_thin++){

          // Generate new proposal and determine it's chi^2
          flag_success=generate_MCMC_chain(MCMC);
          // If this rank is processing chain information ...
          if(my_chain==SID.My_rank){
            // ... and this proposal is kept by the thining then ...
            if(i_thin==n_thin){

              // ... reformat parameter arrays for statistics generation
              for(i_P=0;i_P<n_P;i_P++){
                P_chain_stats[i_P][i_avg%n_avg]=P_chain[i_P];
                P_prop_stats[i_P][i_avg%n_avg] =P_new[i_P];
              }
              ln_Pr_chain[i_avg%n_avg]    =MCMC->ln_Pr_chain;
              ln_Pr_proposals[i_avg%n_avg]=MCMC->ln_Pr_new;

              // ... and write to the chain file
              if(my_chain==SID.My_rank){
                fwrite(&flag_success,             sizeof(char),    1,fp_chain);
                fwrite(&(MCMC->ln_likelihood_new),sizeof(double),  1,fp_chain);
                fwrite(P_new,                     sizeof(double),n_P,fp_chain);
                if(!flag_no_map_write){
                  for(i_DS=0;i_DS<n_DS;i_DS++)
                    fwrite(M_new[i_DS],sizeof(double),n_M[i_DS],fp_chain);
                }
              }
            }
          }

          // Change counters at the end of a thining interval
          if(i_thin==n_thin){
            i_thin=0;
            i_avg++;
          }
        } // i_avg
        if(i_phase==0)
          n_iterations_file_burn++;
        n_iterations_file_total++;
        i_iteration++;

        // Generate statistics for the averaging interval we just completed
        if(my_chain==SID.My_rank){
          // Write the statistics to the chain stats file
          compute_MCMC_chain_stats(P_prop_stats,n_P,n_avg,P_min,P_avg,P_max,dP_avg,auto_cor,slopes,dP_sub,ln_Pr_proposals,&ln_Pr_min,&ln_Pr_avg,&ln_Pr_max);
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
          compute_MCMC_chain_stats(P_chain_stats,n_P,n_avg,P_min,P_avg,P_max,dP_avg,auto_cor,slopes,dP_sub,ln_Pr_chain,&ln_Pr_min,&ln_Pr_avg,&ln_Pr_max);
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
        SID_Bcast(&flag_stop,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
        
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
      SID_log("Proposal success: %5.3f%% (%lld of %lld)",SID_LOG_COMMENT,1e2*((float)(MCMC->n_success)/(float)(MCMC->n_propositions)),MCMC->n_success,MCMC->n_propositions);

      i_iteration      =0;
      i_iteration_start=i_iteration;
      SID_log("Done.",SID_LOG_CLOSE);
    }
    if(my_chain==SID.My_rank){
      fp_chain_config=fopen(filename_chain_config,"wb");
      fwrite(&n_iterations_file_total,sizeof(int),   1,      fp_chain_config);
      fwrite(&n_iterations_file_burn, sizeof(int),   1,      fp_chain_config);
      fwrite(&(MCMC->temperature),    sizeof(double),1,      fp_chain_config);
      fwrite(MCMC->V,                 sizeof(double),n_P*n_P,fp_chain_config);
      fclose(fp_chain_config);
      fclose(fp_chain);
      fclose(fp_stats);
    }
    SID_log("Done.",SID_LOG_CLOSE);

  // Perform analysis on chain(s)
  if(MCMC->flag_analysis_on)
      analyze_MCMC(MCMC);

  // Clean-up
  SID_log("Cleaning-up...",SID_LOG_OPEN);
  SID_free(SID_FARG P_min);
  SID_free(SID_FARG P_max);
  SID_free(SID_FARG P_avg);
  SID_free(SID_FARG dP_avg);
  SID_free(SID_FARG ln_Pr_chain);
  SID_free(SID_FARG ln_Pr_proposals);
  SID_free(SID_FARG drift);
  SID_free(SID_FARG slopes);
  SID_free(SID_FARG dP_sub);
  for(i_P=0;i_P<n_P;i_P++){
    SID_free(SID_FARG P_chain_stats[i_P]);
    SID_free(SID_FARG P_prop_stats[i_P]);
    if(flag_autocor_on)
      SID_free(SID_FARG auto_cor[i_P]);
  }
  SID_free(SID_FARG P_chain_stats);
  SID_free(SID_FARG P_prop_stats);
  if(flag_autocor_on)
    SID_free(SID_FARG auto_cor);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}
