#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  double a_min;
  double a_max;
  double z_min;
  double z_max;
  int    n_snap;
  z_min =atof(argv[1]);
  z_max =atof(argv[2]);
  n_snap=atoi(argv[3]);
  a_min =a_of_z(z_max);
  a_max =a_of_z(z_min);

  ADaPS *cosmo;
  init_cosmo_std(&cosmo);

  SID_log("Constructing %d snapshot expansion factors between z=%lf->%lf...",SID_LOG_OPEN,n_snap,z_min,z_max);

  int     i_a=0;
  int     n_a;
  double  a_1,a_2,da;
  double *a_list;
  double *t_list;
  double  dt;
  n_a   =200;
  da    =(a_max-a_min)/(double)(n_a-1);
  a_list=(double *)SID_malloc(sizeof(double)*n_a);
  t_list=(double *)SID_malloc(sizeof(double)*n_a);
  a_1   =a_min;
  a_2   =a_1+da;
  a_list[i_a]=a_min;
  t_list[i_a]=0.;
  i_a++;
  for(;i_a<n_a;i_a++){
     dt         =deltat_a(&cosmo,a_2,a_1);
     a_list[i_a]=a_2;
     t_list[i_a]=t_list[i_a-1]+dt;
     a_1 =a_2;
     if(i_a==(n_a-2))
       a_2=a_max;
     else
       a_2+=da;
  }
  for(i_a=0;i_a<n_a;i_a++)
     t_list[i_a]/=(1e9*S_PER_YEAR);

  double delta_t;
  delta_t=t_list[n_a-1];
  SID_log("delta_t=%le Gyrs",SID_LOG_COMMENT,delta_t);

  interp_info *interp;
  init_interpolate(t_list,a_list,n_a,gsl_interp_cspline,&interp);

  double a_snap;
  double t_snap;
  for(i_a=0,t_snap=0.;i_a<n_snap;i_a++,t_snap+=delta_t/(double)(n_snap-1)){
     if(i_a==0)
        a_snap=a_min;
     else if(i_a==(n_snap-1))
        a_snap=a_max;
     else
        a_snap=interpolate(interp,t_snap);
     printf("%lf\n",a_snap);
  }  

  SID_free(SID_FARG a_list);
  SID_free(SID_FARG t_list);
  free_interpolate(SID_FARG interp,NULL);

  SID_log("Done.",SID_LOG_CLOSE);

  return(ERROR_NONE);
}

