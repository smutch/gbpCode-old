#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>

double dn_dyn_dt(double a,void *cosmo_in);
double dn_dyn_dt(double a,void *cosmo_in){
  cosmo_info *cosmo=(cosmo_info *)cosmo_in;
  return(1./(a*H_convert(H_z(z_of_a(a),cosmo))*t_dyn_z(z_of_a(a),cosmo)));
}

double delta_n_dyn(double a_1,double a_2,cosmo_info **cosmo);
double delta_n_dyn(double a_1,double a_2,cosmo_info **cosmo){
  // Initialize integral
  double integral;
  int    n_int       =512;
  double abs_accuracy=0.;
  double rel_accuracy=1e-4;
  double abs_error;
  double a_lo        =MIN(a_1,a_2);
  double a_hi        =MAX(a_1,a_2);
  gsl_integration_workspace *wspace;
  gsl_function               integrand;
  integrand.function=dn_dyn_dt;
  integrand.params  =(void *)(*cosmo);
  wspace            =gsl_integration_workspace_alloc(16*n_int);
  gsl_integration_qag(&integrand,
                      a_lo,a_hi,
                      abs_accuracy,rel_accuracy,
                      n_int,
                      GSL_INTEG_GAUSS61,
                      wspace,
                      &integral,&abs_error);
  gsl_integration_workspace_free(wspace);
  return(integral);
}

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL,NULL);

  // Parse arguments
  char   filename_cosmo[MAX_FILENAME_LENGTH];
  strcpy(filename_cosmo,   argv[1]);
  double z_1 =(double)atof(argv[2]);
  double z_2 =(double)atof(argv[3]);
  int    n_z =(int)   atoi(argv[4]);
  double z_lo=MIN(z_1,z_2);
  double z_hi=MAX(z_1,z_2);

  SID_log("Creating snapshot list...",SID_LOG_OPEN);

  // Initialize cosmology
  cosmo_info *cosmo=NULL;
  read_gbpCosmo_file(&cosmo,filename_cosmo);

  // Creeate interpolation tables
  int     n_table    =500;
  double  a_lo_table =a_of_z(z_hi); 
  double  a_hi_table =a_of_z(z_lo); 
  double *a_table    =(double *)SID_malloc(n_table*sizeof(double));
  double *t_table    =(double *)SID_malloc(n_table*sizeof(double));
  double *n_dyn_table=(double *)SID_malloc(n_table*sizeof(double));
  double  da         =(a_hi_table-a_lo_table)/(double)(n_table-1);
  for(int i_table=0;i_table<n_table;i_table++){
     // Expansion factor
     if(i_table==0)
        a_table[i_table]=a_lo_table;
     else if(i_table==(n_table-1))
        a_table[i_table]=a_of_z(z_lo);
     else 
        a_table[i_table]=a_table[0]+da*(double)i_table;
     // Cosmic time
     t_table[i_table]=t_age_a(a_table[i_table],&cosmo);
     // Number of dynamical times
     if(i_table==0)
        n_dyn_table[i_table]=0.;
     else
        n_dyn_table[i_table]=delta_n_dyn(a_table[0],a_table[i_table],&cosmo);
  }
  interp_info *a_of_n_interp=NULL;
  interp_info *t_of_n_interp=NULL;
  init_interpolate(n_dyn_table,a_table,n_table,gsl_interp_akima,&a_of_n_interp);
  init_interpolate(n_dyn_table,t_table,n_table,gsl_interp_akima,&t_of_n_interp);

  // Generate desired table
  double a_last =0.;
  double z_last =0.;
  double t_last =0.;
  double a_i    =0.;
  double z_i    =0.;
  double t_i    =0.;
  double dt     =0.;
  double dz     =0.;
  double n_dyn  =0.;
  double dn_dyn =(n_dyn_table[n_table-1]-n_dyn_table[0])/(double)(n_z-1);
  int i_column=1;
  printf("# Table of expansion factors spaced equally in dynamical time intervals between z=%.2lf and %.2lf\n",z_lo,z_hi);
  printf("# Column (%02d): a:=expansion factor\n",          i_column++);
  printf("#        (%02d): da\n",                           i_column++);
  printf("#        (%02d): z:=redshift\n",                  i_column++);
  printf("#        (%02d): dz\n",                           i_column++);
  printf("#        (%02d): t:=cosmic time  [years]\n",      i_column++);
  printf("#        (%02d): dynamical time [years]\n",       i_column++);
  printf("#        (%02d): dt\n",                           i_column++);
  printf("#        (%02d): n_dyn:=No. of dynamical times\n",i_column++);
  printf("#        (%02d): dn_dyn\n",                       i_column++);
  for(int i_z=0;i_z<n_z;i_z++,n_dyn+=dn_dyn){
     z_last=z_i;
     a_last=a_i;
     t_last=t_i;
     a_i   =interpolate(a_of_n_interp,n_dyn);
     z_i   =z_of_a(a_i);
     t_i   =interpolate(t_of_n_interp,n_dyn)/S_PER_YEAR;
     da    =a_i-a_last;
     dz    =fabs(z_i-z_last);
     dt    =t_i-t_last;
     printf("%le %le %le %le %le %le %le %le %le\n",a_i,da,z_i,dz,t_i,t_dyn_z(z_of_a(a_i),cosmo)/S_PER_YEAR,dt,n_dyn,dn_dyn);
  }

  // Clean-up
  SID_free(SID_FARG a_table);
  SID_free(SID_FARG t_table);
  SID_free(SID_FARG n_dyn_table);
  free_interpolate(SID_FARG a_of_n_interp,NULL);
  free_interpolate(SID_FARG t_of_n_interp,NULL);
  free_cosmo(&cosmo);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

