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

void read_MCMC_state(MCMC_info *MCMC){
  char      filename_output_dir[MAX_FILENAME_LENGTH];
  char      filename_chain_dir[MAX_FILENAME_LENGTH];
  char      filename_results_dir[MAX_FILENAME_LENGTH];
  char      filename_plots_dir[MAX_FILENAME_LENGTH];
  char      filename_run[MAX_FILENAME_LENGTH];
  char      filename_chain[MAX_FILENAME_LENGTH];
  char      filename_chain_config[MAX_FILENAME_LENGTH];
  char      filename_stats[MAX_FILENAME_LENGTH];
  char      filename_coverage[MAX_FILENAME_LENGTH];
  char      filename_chain_covariance[MAX_FILENAME_LENGTH];
  char      filename_covariance[MAX_FILENAME_LENGTH];
  char      filename_histograms[MAX_FILENAME_LENGTH];
  char      filename_results[MAX_FILENAME_LENGTH];
  char      filename_stop[MAX_FILENAME_LENGTH];
  char      format_string[32];
  int       my_chain;
  int       i_P,i_DS,i_M,i_array;
  double   *V_read;
  FILE     *fp_run;
  FILE     *fp_chain;
  FILE     *fp_chain_config;
  FILE     *fp_stats;
  FILE     *fp_coverage;
  FILE     *fp_chain_covariance;
  FILE     *fp_covariance;
  FILE     *fp_histograms;
  FILE     *fp_results;
  FILE     *fp_stop;
  MCMC_DS_info *current_DS;

  set_MCMC_mode(MCMC,MCMC_MODE_DEFAULT);
  my_chain=MCMC->my_chain;

    SID_log("Reading MCMC state from {%s}...",SID_LOG_OPEN,MCMC->filename_output_dir);

    // Set directories
    sprintf(filename_output_dir, "%s/",        MCMC->filename_output_dir);
    sprintf(filename_chain_dir,  "%s/chains/", MCMC->filename_output_dir);
    sprintf(filename_results_dir,"%s/results/",MCMC->filename_output_dir);
    sprintf(filename_plots_dir,  "%s/plots/",  MCMC->filename_output_dir);
    // Set filenames
    sprintf(filename_run,             "%s/run.dat",                  MCMC->filename_output_dir);
    sprintf(filename_chain,           "%s/chain_trace_%06d.dat",     filename_chain_dir,my_chain);
    sprintf(filename_chain_config,    "%s/chain_config_%06d.dat",    filename_chain_dir,my_chain);
    sprintf(filename_chain_covariance,"%s/chain_covariance_%06d.dat",filename_chain_dir,my_chain);
    sprintf(filename_stats,           "%s/chain_stats_%06d.dat",     filename_chain_dir,my_chain);
    sprintf(filename_coverage,        "%s/coverage.dat",             filename_results_dir);
    sprintf(filename_histograms,      "%s/histograms.dat",           filename_results_dir);
    sprintf(filename_covariance,      "%s/covariance.dat",           filename_results_dir);

    MCMC->map_P_to_M                =NULL;
    MCMC->compute_MCMC_ln_likelihood=compute_MCMC_ln_likelihood_default;
    MCMC->params                    =NULL;
    MCMC->temperature               =1.0;
    MCMC->n_P                       =0;
    MCMC->n_thin                    =1;
    MCMC->n_DS                      =0;
    MCMC->n_M_total                 =0;
    MCMC->n_arrays                  =0;
    MCMC->n_M                       =NULL;
    MCMC->array                     =NULL;
    MCMC->V                         =NULL;
    MCMC->m                         =NULL;
    MCMC->b                         =NULL;
    MCMC->RNG                       =NULL;
    MCMC->flag_integrate_on         =TRUE;
    MCMC->flag_analysis_on          =TRUE;
    MCMC->first_map_call            =TRUE;
    MCMC->mode                      =MCMC_MODE_DEFAULT;
    MCMC->DS                        =NULL;
    MCMC->last                      =NULL;
    MCMC->problem_name              =(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);

    // Read/Write Header file
    if((fp_run=fopen(filename_run,"rb"))!=NULL){
      fp_run=fopen(filename_run,"rb");
      fread(MCMC->problem_name,        sizeof(char),MCMC_NAME_SIZE,fp_run);
      fread(&(MCMC->n_avg),            sizeof(int),              1,fp_run);
      fread(&(MCMC->flag_autocor_on),  sizeof(int),              1,fp_run);
      fread(&(MCMC->flag_no_map_write),sizeof(int),              1,fp_run);
      fread(&(MCMC->n_P),              sizeof(int),              1,fp_run);
      SID_log("Problem name    ={%s}",SID_LOG_COMMENT,MCMC->problem_name);
      SID_log("n_avg           ={%d}",SID_LOG_COMMENT,MCMC->n_avg);
      SID_log("flag_autocor_on ={%d}",SID_LOG_COMMENT,MCMC->flag_autocor_on);
      MCMC->P_names    =(char  **)SID_malloc(sizeof(char *)*MCMC->n_P);
      MCMC->P_init     =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      MCMC->P_new        =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      MCMC->P_last       =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      MCMC->P_chain      =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      MCMC->P_limit_min  =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      MCMC->P_limit_max  =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
      for(i_P=0;i_P<MCMC->n_P;i_P++)
        MCMC->P_limit_min[i_P]=-DBL_MAX*1e-3;
      for(i_P=0;i_P<MCMC->n_P;i_P++)
        MCMC->P_limit_max[i_P]=DBL_MAX*1e-3;

      SID_log("Parameters (name, initial_value,limit min,limit max):",SID_LOG_OPEN);
      MCMC->P_name_length=0;
      for(i_P=0;i_P<MCMC->n_P;i_P++){
        MCMC->P_names[i_P]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
        fread(MCMC->P_names[i_P],       sizeof(char),  MCMC_NAME_SIZE,fp_run);
        fread(&(MCMC->P_init[i_P]),     sizeof(double),             1,fp_run);
        fread(&(MCMC->P_limit_min[i_P]),sizeof(double),             1,fp_run);
        fread(&(MCMC->P_limit_max[i_P]),sizeof(double),             1,fp_run);
        MCMC->P_name_length=MAX(MCMC->P_name_length,strlen(MCMC->P_names[i_P]));
      }
      sprintf(MCMC->P_name_format,"%%-%ds",            MCMC->P_name_length);
      sprintf(format_string,      "%s %%13.6le %%13.6le %%13.6le",MCMC->P_name_format);
      for(i_P=0;i_P<MCMC->n_P;i_P++)
        SID_log(format_string,SID_LOG_COMMENT,MCMC->P_names[i_P],MCMC->P_init[i_P],MCMC->P_limit_min[i_P],MCMC->P_limit_max[i_P]);
      SID_log(NULL,SID_LOG_CLOSE|SID_LOG_NOPRINT);
      fread(&(MCMC->n_arrays),sizeof(int),1,fp_run);
      SID_log("n_arrays=%d",  SID_LOG_OPEN,MCMC->n_arrays);
      MCMC->array     =(double **)SID_malloc(sizeof(double *)*MCMC->n_arrays);
      MCMC->array_name=(char   **)SID_malloc(sizeof(char   *)*MCMC->n_arrays);
      for(i_array=0;i_array<MCMC->n_arrays;i_array++){
        MCMC->array[i_array]     =(double *)SID_malloc(sizeof(double)*MCMC->n_P);
        MCMC->array_name[i_array]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
        fread(MCMC->array_name[i_array],sizeof(char),  MCMC_NAME_SIZE,fp_run);
        fread(MCMC->array[i_array],     sizeof(double),MCMC->n_P,     fp_run);
        SID_log("array #%03d name ={%s}",SID_LOG_COMMENT,i_array,MCMC->array_name[i_array]);
      }
      SID_log(NULL,SID_LOG_CLOSE|SID_LOG_NOPRINT);
      fread(&(MCMC->n_DS),sizeof(int),1,fp_run);
      SID_log("Reading %d datasets...",SID_LOG_OPEN,MCMC->n_DS);
      for(i_DS=0;i_DS<MCMC->n_DS;i_DS++){
        SID_log("Dataset #%03d:",SID_LOG_OPEN,i_DS);
        current_DS           =(MCMC_DS_info *)SID_malloc(sizeof(MCMC_DS_info));
        fread(current_DS->name,  sizeof(char),MCMC_NAME_SIZE,fp_run);
        fread(&(current_DS->n_M),sizeof(int),              1,fp_run);
        MCMC->n_M_total+=current_DS->n_M;
        SID_log("name    ={%s}",SID_LOG_COMMENT,current_DS->name);
        SID_log("n_M     =%d",  SID_LOG_COMMENT,current_DS->n_M);
        current_DS->M_target =(double *)SID_malloc(sizeof(double)*current_DS->n_M);
        current_DS->dM_target=(double *)SID_malloc(sizeof(double)*current_DS->n_M);
        current_DS->params   =NULL;
        fread(current_DS->M_target,   sizeof(double),current_DS->n_M,fp_run);
        fread(current_DS->dM_target,  sizeof(double),current_DS->n_M,fp_run);
        fread(&(current_DS->n_arrays),sizeof(int),                 1,fp_run);
        SID_log("n_arrays=%d",  SID_LOG_OPEN,current_DS->n_arrays);
        current_DS->array     =(double **)SID_malloc(sizeof(double *)*current_DS->n_arrays);
        current_DS->array_name=(char   **)SID_malloc(sizeof(char   *)*current_DS->n_arrays);
        for(i_array=0;i_array<current_DS->n_arrays;i_array++){
          current_DS->array[i_array]     =(double *)SID_malloc(sizeof(double)*current_DS->n_M);
          current_DS->array_name[i_array]=(char *)SID_malloc(sizeof(char)*MCMC_NAME_SIZE);
          fread(current_DS->array_name[i_array],sizeof(char),  MCMC_NAME_SIZE, fp_run);
          fread(current_DS->array[i_array],     sizeof(double),current_DS->n_M,fp_run);
          SID_log("array #%03d name={%s}",SID_LOG_COMMENT,i_array,current_DS->array_name[i_array]);
        }
        SID_log(NULL,SID_LOG_CLOSE|SID_LOG_NOPRINT);
        current_DS->next=NULL;
        if(MCMC->DS==NULL)
          MCMC->DS=current_DS;
        else
          MCMC->last->next=current_DS;
        MCMC->last=current_DS;
        SID_log(NULL,SID_LOG_CLOSE|SID_LOG_NOPRINT);
      }
      SID_log("Done.",SID_LOG_CLOSE);
      fclose(fp_run);
    }

    // ... fetch the number of intervals that have already been computed ...
    fp_chain_config=fopen(filename_chain_config,"rb");
    V_read=(double *)SID_malloc(sizeof(double)*MCMC->n_P*MCMC->n_P);
    fread(&(MCMC->n_iterations),     sizeof(int),   1,      fp_chain_config);
    fread(&(MCMC->n_iterations_burn),sizeof(int),   1,      fp_chain_config);

    // ... fetch the temperature and covariance matrix that was being used
    fread(&(MCMC->temperature),      sizeof(double),1,                  fp_chain_config);
    fread(V_read,                    sizeof(double),MCMC->n_P*MCMC->n_P,fp_chain_config);
    set_MCMC_covariance(MCMC,V_read);
    SID_free(SID_FARG V_read);

    // Initialize dataset arrays
    init_MCMC_arrays(MCMC);

    SID_log("# burn  iterations = %d", SID_LOG_COMMENT,MCMC->n_iterations_burn);
    SID_log("# total iterations = %d", SID_LOG_COMMENT,MCMC->n_iterations);
    SID_log("Temperature        = %le",SID_LOG_COMMENT,MCMC->temperature);
    fclose(fp_chain_config);

    SID_log("Done.",SID_LOG_CLOSE);
}
