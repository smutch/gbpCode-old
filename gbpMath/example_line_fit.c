#define  _MAIN  // This always needs to be at the very top of a file if it has a main() function in it
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>  // Always needed by gbpMCMC
#include <gbpMath.h> // Always needed by gbpMCMC

// This structure holds any information that you want to pass to the model function 
typedef struct model_params model_params;
struct model_params{
  int    param1;
  double param2;
};

// This function must generate and return the model
int model_function(double *P,MCMC_info *MCMC,double **M);
int model_function(double *P,MCMC_info *MCMC,double **M){
  MCMC_DS_info  *current_DS;
  MCMC_DS_info  *next_DS;
  int            i_P,i_DS,i_M;
  double        *x;
  static int     param1;
  static double  param2;

  // These variables get set at the first call and remain set
  //   Although param1 and param2 do not get used in this
  //   example, it shows how to pass information from the
  //   main function to the model function
  if(MCMC->first_map_call){
    param1=((model_params *)(MCMC->params))->param1;
    param2=((model_params *)(MCMC->params))->param2;
  }
  // param1 and param2 are set for all calls of the
  //   model function now.  n.b.: they MUST be
  //   declared as static variables for this to work

  // COMPUTE THE MODEL

  // Loop over each dataset ...
  i_DS=0;
  i_P =0;
  current_DS=MCMC->DS;
  while(current_DS!=NULL){
    next_DS=current_DS->next;
    // You can store as many arrays for each dataset as you want.  In this
    //   case, we have stored the x-values of the dataset in the first array.
    x=current_DS->array[0];
    // Loop over each element of the dataset ...
    for(i_M=0;i_M<current_DS->n_M;i_M++)
      M[i_DS][i_M]=P[i_P+0]*x[i_M]+P[i_P+1]; // Even-index parameters are slopes and
                                             //   odd-indexed parameters are intercepts
                                             //   for this example
    current_DS=next_DS;
    i_DS++;
    i_P+=2;
  }

  // If the model is invalid for some reason,
  //   return MCMC_MAP_RETURN_BAD instead and
  //   a new parameter set will be chosen and
  //   this model function reacalled until
  //   MCMC_MAP_RETURN_GOOD is returned.
  return(MCMC_MAP_RETURN_GOOD);
}

int main(int argc, char *argv[]){
  int          i_DS;
  int          i_M;
  int          n_DS;
  int          i_P;
  int          n_M;
  model_params params;
  MCMC_info    MCMC;
  RNG_info     RNG;

  // This must be called at the very beginning of the main function  
  SID_init(&argc,&argv,NULL,NULL);

  SID_log("Running the MCMC example code '%s'...",SID_LOG_OPEN,argv[0]);

  // Parse the parameters passed to the code from stdin.  On this case, 
  //   read the number of datasets and the dataset size 
  SID_log("Parsing the input parameters...",SID_LOG_OPEN);
  if(argc==3){
    n_DS=atoi(argv[1]);
    n_M =atoi(argv[2]);
    SID_log("%d datasets of size %d specified...",SID_LOG_CONTINUE,n_DS,n_M);
  }
  else{
    SID_log("SYNTAX: example_fit_line n_datasets size_of_datasets",SID_LOG_COMMENT);
    SID_exit(ERROR_SYNTAX);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize MCMC params structure (not really used for this example)
  //   This structure allows you to pass information to the model function if need-be.
  params.param1=0;
  params.param2=0.;

  // Initialize a random number generator for creating the example datasets
  int seed;
  seed=-1; // Initialize RNG from the clock
  init_RNG(&seed,&RNG,RNG_DEFAULT);

  // Create an array of strings which holds the names of the parameters
  char   **P_names;
  double  *P_init;
  double  *P_limit_min;
  double  *P_limit_max;
  SID_log("Creating the actual values the fit should give us...",SID_LOG_OPEN);
  P_names    =(char   **)SID_malloc(sizeof(char   *)*2*n_DS);
  P_init     =(double  *)SID_malloc(sizeof(double *)*2*n_DS);
  P_limit_min=(double  *)SID_malloc(sizeof(double *)*2*n_DS);
  P_limit_max=(double  *)SID_malloc(sizeof(double *)*2*n_DS);
  for(i_DS=0,i_P=0;i_DS<n_DS;i_DS++,i_P+=2){
    SID_log("Dataset #%02d:",SID_LOG_OPEN,i_DS+1);
    // Initialize slope of random example dataset
    P_names[i_P+0]    =(char *)SID_malloc(sizeof(char)*32);
    sprintf(P_names[i_P+0],"slope_%d",i_DS);
    P_init[i_P+0]     =random_number(&RNG);
    P_limit_min[i_P+0]=0.;
    P_limit_max[i_P+0]=1.;
    SID_log("Slope    =%f",SID_LOG_COMMENT,P_init[i_P+0]);
    // Initialize intercept of random example dataset
    P_names[i_P+1]    =(char *)SID_malloc(sizeof(char)*32);
    sprintf(P_names[i_P+1],"intercept_%d",i_DS);
    P_init[i_P+1]     =random_number(&RNG);
    P_limit_min[i_P+1]=0.;
    P_limit_max[i_P+1]=1.;
    SID_log("Intercept=%f",SID_LOG_COMMENT,P_init[i_P+1]);
    SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize the MCMC calculation
  init_MCMC(&MCMC,
            "This is the example_line_fit MCMC Run Title",
            (void *)(&params),   // This is a void pointer that can be used to pass any kind of information to the model
            &model_function,     // This is the name of the function that will be called to generate the model
            2*n_DS,              // The number of parameters in the fit.  We'll have a slope and an intercept per data set.
            P_init,              // These will be the starting values for the fit. (Also the correct values for this example)
            P_names,             // This array of strings holds the parameter names
            P_limit_min,         // This array holds the minimum values allowed by the fit
            P_limit_max,         // This array holds the maximum values allowed by the fit
            0);                  // Optionally we could pass through some arrays here if we wished

  // Loop over the number of datasets
  double *ordinate_values;
  double *data_values;
  double *uncertainties;
  char    dataset_name[32];
  ordinate_values=(double *)SID_malloc(sizeof(double)*n_M);
  data_values    =(double *)SID_malloc(sizeof(double)*n_M);
  uncertainties  =(double *)SID_malloc(sizeof(double)*n_M);
  for(i_DS=0,i_P=0;i_DS<n_DS;i_DS++,i_P+=2){
    // Set the dataset name
    sprintf(dataset_name,"Dataset #%02d",i_DS+1);
    SID_log("Generating %s and adding it to the calculation...",SID_LOG_OPEN,dataset_name);

    // Create the i_DS'th data set
    for(i_M=0;i_M<n_M;i_M++){
      ordinate_values[i_M]=(double)random_number(&RNG);                      // Return a value from [0->1)
      data_values[i_M]    =P_init[i_P+0]*ordinate_values[i_M]+P_init[i_P+1]; // Initially put the data value on the line
      // Create a positive uncertainty.  If an uncertainty is <=0,
      //   the MCMC fit will get caught in an infinite loop
      uncertainties[i_M]=0.;
      while(uncertainties[i_M]<=0.)
        uncertainties[i_M]=data_values[i_M]*(double)random_number(&RNG)/5.;     // Choose an uncertainty from (0% to 20%)
      data_values[i_M]+=uncertainties[i_M]*(1.-2.*(double)random_number(&RNG)); // Scatter the data value by this uncertainty
      SID_log("%10.3lf %10.3lf %10.3lf",SID_LOG_COMMENT,ordinate_values[i_M],data_values[i_M],uncertainties[i_M]);
    }

    // Add the i_DS'th data set to the fit.  The arrays get copied in add_MCMC_DS() so we can reuse them here.
    add_MCMC_DS(&MCMC,
                dataset_name,           // The name of the dataset
                1,                      // The dimensionality of the dataset (n_D)
                &n_M,                   // An n_D-long array giving the dimension sizes of the dataset
                data_values,            // The data values of the dataset
                uncertainties,          // The uncertainties of the dataset
                NULL,                   // A pointer to pass information to the dataset
                1,                      // The number of arrays associated with this dataset
                ordinate_values,        // Values for the first array: An array holding the ordinate values of this dataset
                "x-values [unitless]"); // Name   for the first array

    SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up arrays (not needed anymore because they get copied to the MCMC structure)
  for(i_DS=0,i_P=0;i_DS<n_DS;i_DS++){
     SID_free(SID_FARG P_names[i_P++]);
     SID_free(SID_FARG P_names[i_P++]);
  }
  SID_free(SID_FARG P_names);
  SID_free(SID_FARG P_init);
  SID_free(SID_FARG P_limit_min);
  SID_free(SID_FARG P_limit_max);
  SID_free(SID_FARG ordinate_values);
  SID_free(SID_FARG data_values);
  SID_free(SID_FARG uncertainties);

  // Free the random number generator used to generate the example datasets
  free_RNG(&RNG); 
  
  // This sets the directory that will be created to store the calculation
  set_MCMC_directory(&MCMC,"example_line_fit_MCMC");

  // Set the parameters that determine the performance and accuracy of the calculation.
  int    n_avg;
  int    n_thin;
  int    n_burn;
  int    n_integrate;
  int    coverage_size;
  double success_target;
  double success_threshold;
  double covariance_threshold;
  int    n_autotune;
  int    n_autotune_randomize;
  int    n_autotune_temperature;
  int    n_autotune_covariance;
  n_avg           =        1000; // Chains are calculated in "averaging intervals" (ie. chunks) of this size.
                                 //  Should be at least 100.
  n_thin          =           1; // Every n_thin'th element of the chain is discarded.
                                 //  Usually set to 1.
  n_burn          =         100; // This is the number of averaging intervals that are discarded for the burn-in.
  n_integrate     =        1000; // This is the number of averaging intervals that are used for integration.
                                 //  The final chains will be of size n_avg*n_integrate
  coverage_size   =         100; // This is the number of bins/pixel-dimension used to create histograms/degeneracy plots.
  success_target  =         35.; // This is the desired success rate for propositions (in percent).
                                 //  Should be about 35%
  success_threshold=         5.; // When performing autotuning, this is the threshold within-which we demand
                                 //  the success rate to converge (in percent).
                                 //  Should be about 5%.
  covariance_threshold=     50.; // When performing autotuning, this is the threshold within-which we demand
                                 //  the covariance matrix to converge (in percent).
                                 //  Ideally, this is ~10% but if there is very little covariance
                                 //  between parameters, it may need to be made very large (eg. 1e7)
  n_autotune          =       2; // The number of passes to be performed for the autotuning.
                                 //  Should be at least 2.  Rarely needs to be more than 3.
  n_autotune_randomize=       0; // The number of model calls that will be used to randomize the starting position.
                                 //  Really only needs to be >0 if you are running multiple chains.
  n_autotune_temperature= 10000; // The number of model calls that will be used to tune the chain temperature.
                                 //  Should be AT LEAST 100 ... ideally 1000.
  n_autotune_covariance=1000000; // The number of model calls that will be used to tune the covariance matrix.
                                 //  Ideally this is pretty large, especially if there is little covariance
                                 //  between some parameters.  Otherwise, tuning may never converge.

  // Set the integration parameters of the calculation
  set_MCMC_iterations(&MCMC,
                      n_avg,
                      n_thin,
                      n_burn,
                      n_integrate);
  set_MCMC_autotune(&MCMC,
                    success_target,
                    success_threshold,
                    covariance_threshold,
                    n_autotune,
                    n_autotune_randomize,
                    n_autotune_temperature,
                    n_autotune_covariance);

  // Perform the MCMC calculation
  compute_MCMC(&MCMC);

  // At this point, the MCMC calculation is done and the
  //   results of the fit will be in the 'results' directory
  //   of the MCMC output directory.

  // Clean-up
  free_MCMC(&MCMC);

  // Exit
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);

  // Once the code exits, run the 'allresults_MCMC.py'
  //   script on the directory and it will generate 
  //   all the diagnostic plots.  You need to have 
  //   python installed for it to work though.
}

