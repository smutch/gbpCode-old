#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_linear_theory.h>
#include <gsl/gsl_multimin.h>

typedef struct init_gbpCosmo2gbpCosmo_integrand_params_struct init_gbpCosmo2gbpCosmo_integrand_params;
struct init_gbpCosmo2gbpCosmo_integrand_params_struct {
  double       inv_s;
  double       z_source;
  double       z_target;
  double       R_1;
  double       R_2;
  cosmo_info **cosmo_source;
  cosmo_info **cosmo_target;
  int                        n_int;
  gsl_integration_workspace *wspace;
};

double init_gbpCosmo2gbpCosmo_integrand(double R,void *params_in){
  init_gbpCosmo2gbpCosmo_integrand_params *params =(init_gbpCosmo2gbpCosmo_integrand_params *)params_in;
  double       inv_s       =params->inv_s;
  double       z_source      =params->z_source;
  double       z_target    =params->z_target;
  cosmo_info **cosmo_source  =params->cosmo_source;
  cosmo_info **cosmo_target=params->cosmo_target;
  return(pow(1.-(sigma_R(cosmo_source,inv_s*R,z_source,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER))/(sigma_R(cosmo_target,R,z_target,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER)),2.)/R);
}

double init_gbpCosmo2gbpCosmo_minimize_function(const gsl_vector *v_i,void *params_in){
   // Set variable integrand parameters
   init_gbpCosmo2gbpCosmo_integrand_params *params=(init_gbpCosmo2gbpCosmo_integrand_params *)params_in;
   params->inv_s   =gsl_vector_get(v_i,0);
   params->z_target=gsl_vector_get(v_i,1);

   // Perform integral to minimize
   gsl_function integrand;
   double       delta_i;
   double       abs_error;
   integrand.function=init_gbpCosmo2gbpCosmo_integrand;
   integrand.params  =params_in;
   gsl_integration_qag(&integrand,
                       params->R_1,params->R_2,
                       0,1e-3,
                       params->n_int,
                       GSL_INTEG_GAUSS61,
                       params->wspace,
                       &delta_i,&abs_error);
   return(delta_i/take_ln(params->R_2/params->R_1));
}

void init_gbpCosmo2gbpCosmo(cosmo_info      **cosmo_source,
                            cosmo_info      **cosmo_target,
                            double            z_min,
                            double            M_min,
                            double            M_max,
                            gbpCosmo2gbpCosmo_info *gbpCosmo2gbpCosmo){
   SID_log("Initializing cosmology scaling...",SID_LOG_OPEN|SID_LOG_TIMER);
   SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

   // Store some infor in the gbpCosmo2gbpCosmo_info structure
   gbpCosmo2gbpCosmo->M_min       =M_min;
   gbpCosmo2gbpCosmo->M_max       =M_max;
   gbpCosmo2gbpCosmo->z_min       =z_min;
   gbpCosmo2gbpCosmo->cosmo_source=(*cosmo_source);
   gbpCosmo2gbpCosmo->cosmo_target=(*cosmo_target);

   // Perform minimization
   //const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
   const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex;
   gsl_multimin_fminimizer            *s = NULL;
   gsl_vector *ss, *x;
   gsl_multimin_function minex_func;
 
   // Starting point 
   x = gsl_vector_alloc (2);
   gsl_vector_set (x, 0, 1.);    // inv_s
   gsl_vector_set (x, 1, z_min); // z_scaled
 
   // Set initial step sizes to 1 
   ss = gsl_vector_alloc (2);
   gsl_vector_set_all (ss, 1.0);

   // Set parameters
   init_gbpCosmo2gbpCosmo_integrand_params params;
   params.cosmo_source=cosmo_source;
   params.cosmo_target=cosmo_target;
   params.z_source    =z_min;
   params.R_1         =R_of_M(M_min,*cosmo_source);
   params.R_2         =R_of_M(M_max,*cosmo_source);
   params.inv_s       =gsl_vector_get(x,0);
   params.z_target    =gsl_vector_get(x,1);
   params.n_int       =100;
   params.wspace      =gsl_integration_workspace_alloc(params.n_int);

   // Initialize method
   minex_func.n      = 2;
   minex_func.f      = init_gbpCosmo2gbpCosmo_minimize_function;
   minex_func.params = (void *)(&params);
   s                 = gsl_multimin_fminimizer_alloc (T, 2);
   gsl_multimin_fminimizer_set(s,&minex_func,x,ss);

   // Perform minimization 
   double size;
   int    status;
   size_t iter    =  0;
   size_t iter_max=200;
   do{
       iter++;
       status=gsl_multimin_fminimizer_iterate(s);
       if(status) 
          SID_trap_error("Error encountered during minimisation in init_gbpCosmo2gbpCosmo() (status=%d).",ERROR_LOGIC,status);
       size   = gsl_multimin_fminimizer_size(s);
       status = gsl_multimin_test_size(size,1e-2);
   } while(status==GSL_CONTINUE && iter<=iter_max);
   if(status!=GSL_SUCCESS)
      SID_trap_error("Failed to converge during minimisation in init_gbpCosmo2gbpCosmo() (status=%d,iter=%d).",ERROR_LOGIC,status,iter);

   // Finalize results   
   double Omega_M_source =    ((double *)ADaPS_fetch(*cosmo_source,"Omega_M") )[0];
   double H_Hubble_source=1e2*((double *)ADaPS_fetch(*cosmo_source,"h_Hubble"))[0];
   double Omega_M_target =    ((double *)ADaPS_fetch(*cosmo_target,"Omega_M") )[0];
   double H_Hubble_target=1e2*((double *)ADaPS_fetch(*cosmo_target,"h_Hubble"))[0];
   gbpCosmo2gbpCosmo->s_L         =1./gsl_vector_get(s->x,0);
   gbpCosmo2gbpCosmo->s_M         =(Omega_M_target*H_Hubble_target)/(Omega_M_source*H_Hubble_source)*pow((gbpCosmo2gbpCosmo->s_L),3.);
   gbpCosmo2gbpCosmo->z_min_scaled=gsl_vector_get(s->x,1);;

   // Calculate growth factors needed for
   //    determining redshift mappings
   gbpCosmo2gbpCosmo->D_prime_z_min=linear_growth_factor(z_min,                    *cosmo_target);
   gbpCosmo2gbpCosmo->D_z_scaled   =linear_growth_factor(gbpCosmo2gbpCosmo->z_min_scaled,*cosmo_source);
   gbpCosmo2gbpCosmo->D_ratio      =gbpCosmo2gbpCosmo->D_prime_z_min/gbpCosmo2gbpCosmo->D_z_scaled;

   // Clean-up
   gsl_vector_free(x);
   gsl_vector_free(ss);
   gsl_multimin_fminimizer_free(s);
   gsl_integration_workspace_free(params.wspace);
   SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
   SID_log("Done.",SID_LOG_CLOSE);
}
