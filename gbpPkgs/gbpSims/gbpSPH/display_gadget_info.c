#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>

void display_gadget_info(plist_info  *plist){
  FILE   *fp;
  char  **pname;
  size_t  i,j,k; 
  int     counter;
  size_t  n_of_type[N_GADGET_TYPE];
  unsigned int n_of_type_tmp[N_GADGET_TYPE];
  int     flag_used[N_GADGET_TYPE];
  size_t  n_particles;
  double  mass_array[N_GADGET_TYPE];
  int     unused[256];
  size_t  n_all[N_GADGET_TYPE];
  int     n_all_tmp[N_GADGET_TYPE];
  int     n_files;
  int     junk;
  double  d_value;
  double *d_array;
  double *d1_array;
  double *d2_array;
  double *d3_array;
  size_t  min_i;
  size_t  max_i;
  double  mean;
  double  min;
  double  max;
  double  median;
  double  std_dev;
  float   f_temp;
  float   f1_temp;
  float   f2_temp;
  float   f3_temp;
  int     i1_temp;
  int     i2_temp;
  char    var_name[256];
  char    var1_name[256];
  char    var2_name[256];
  char    var3_name[256];
  int     n_type_used;
  long    record_length;
  int     n_return;
  int     s_load;
  int     flag_alloc_d1_array;
  int     i_rank;

  double  redshift;
  double  h_Hubble;
  double  Omega_M;
  double  Omega_Lambda;
  double  rho_crit;

  pname=plist->species;

  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
    if(i_rank==SID.My_rank){
      if(SID.n_proc>1)
        fprintf(stderr,"\nInfo for rank=%d:\n--------------\n",SID.My_rank);

      // Display header
      fprintf(stderr,"\n");
      fprintf(stderr,"Using length   unit=%le kpc\n",  plist->length_unit/M_PER_KPC);
      fprintf(stderr,"      mass     unit=%le M_SOL\n",plist->mass_unit/M_SOL);
      fprintf(stderr,"      velocity unit=%le km/s\n", plist->velocity_unit/1e3);
      fprintf(stderr,"\n");
      fprintf(stderr,"Header:\n------\n");
      display_gadget_header(plist);

      // Determine particle numbers for each species
      n_particles=0;
      for(i=0,n_type_used=0;i<plist->n_species;i++) {
        if(ADaPS_exist(plist->data,"n_%s",pname[i])){
          n_of_type[i]=((size_t *)ADaPS_fetch(plist->data,"n_%s",pname[i]))[0];
          if(n_of_type[i]>0){
            n_particles+=n_of_type[i];
            flag_used[i]=TRUE;
            n_type_used++;
          }
          else{
            n_of_type[i]=0;
            flag_used[i]=FALSE;
          }
        }
        else{
          n_of_type[i]=0;
          flag_used[i]=FALSE;
        }
      }

      /* Don't write anything more than the header if there's no particles */
      if(n_particles>0){

        h_Hubble=((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];

        fprintf(stderr,"\n");
        fprintf(stderr,"%11s  %11s  %11s  %11s  %11s  %11s\n",  "Variable","Min","Max","Mean","Median","Std. Dev'n");
        fprintf(stderr,"%11s  %11s  %11s  %11s  %11s  %11s\n\n","--------","---","---","----","------","----------");

        // Positions
        for(i=0;i<plist->n_species;i++) {
          if(flag_used[i]){
            sprintf(var_name,"x_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=(M_PER_KPC/h_Hubble);
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=(M_PER_KPC/h_Hubble);
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=(M_PER_KPC/h_Hubble);
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=(M_PER_KPC/h_Hubble);
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=(M_PER_KPC/h_Hubble);
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le kpc/h\n",var_name,min,max,mean,median,std_dev);
            }

            sprintf(var_name,"y_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=(M_PER_KPC/h_Hubble);
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=(M_PER_KPC/h_Hubble);
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=(M_PER_KPC/h_Hubble);
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=(M_PER_KPC/h_Hubble);
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=(M_PER_KPC/h_Hubble);
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le kpc/h\n",var_name,min,max,mean,median,std_dev);
            }
        
            sprintf(var_name,"z_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=(M_PER_KPC/h_Hubble);
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=(M_PER_KPC/h_Hubble);
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=(M_PER_KPC/h_Hubble);
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=(M_PER_KPC/h_Hubble);
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=(M_PER_KPC/h_Hubble);
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le kpc/h\n",var_name,min,max,mean,median,std_dev);
            }
          }
        }

        // Velocities 
        for(i=0;i<plist->n_species;i++) {
          if(flag_used[i]){
            sprintf(var_name,"vx_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=1e3;
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=1e3;
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=1e3;
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=1e3;
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=1e3;
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le km/s\n",var_name,min,max,mean,median,std_dev);
            }

            sprintf(var_name,"vy_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=1e3;
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=1e3;
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=1e3;
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=1e3;
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=1e3;
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le km/s\n",var_name,min,max,mean,median,std_dev);
            }

            sprintf(var_name,"vz_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=1e3;
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=1e3;
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=1e3;
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=1e3;
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=1e3;
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le km/s\n",var_name,min,max,mean,median,std_dev);
            }

          }
        }

        // Masses 
        for(i=0;i<plist->n_species;i++) {
          if(flag_used[i]){
            sprintf(var_name,"M_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_DOUBLE,CALC_MODE_RETURN_DOUBLE);mean   /=M_SOL;
              calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_DOUBLE,CALC_MODE_RETURN_DOUBLE);median /=M_SOL;
              calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_DOUBLE,CALC_MODE_RETURN_DOUBLE);min    /=M_SOL;
              calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_DOUBLE,CALC_MODE_RETURN_DOUBLE);max    /=M_SOL;
              calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_DOUBLE,CALC_MODE_RETURN_DOUBLE);std_dev/=M_SOL;
              fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le M_sol\n",var_name,min,max,mean,median,std_dev);
            }
          }
        }

        /* Gas properties */
        if(n_of_type[GADGET_TYPE_GAS]>0){
          if(ADaPS_exist(plist->data,"redshift"))
            redshift=((double *)ADaPS_fetch(plist->data,"redshift"))[0];
          else
            redshift=0.;
          h_Hubble    =((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];
          Omega_M     =((double *)ADaPS_fetch(plist->data,"Omega_M"))[0];
          Omega_Lambda=((double *)ADaPS_fetch(plist->data,"Omega_Lambda"))[0];
          rho_crit=rho_crit_z_strip(redshift,h_Hubble,Omega_M,Omega_Lambda);
          sprintf(var_name,"u_%s",pname[GADGET_TYPE_GAS]);
          if(ADaPS_exist(plist->data,var_name)){
            calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=1e3;
            calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=1e3;
            calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=1e3;
            calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=1e3;
            calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=1e3;
            fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le Joules/kg\n",var_name,min,max,mean,median,std_dev);
          }
          sprintf(var_name,"rho_%s",pname[GADGET_TYPE_GAS]);
          if(ADaPS_exist(plist->data,var_name)){
            calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);mean   /=rho_crit;
            calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);median /=rho_crit;
            calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);min    /=rho_crit;
            calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);max    /=rho_crit;
            calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);std_dev/=rho_crit;
            fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le [rho_crit(z)]\n",var_name,min,max,mean,median,std_dev);
          }
          sprintf(var_name,"T_%s",pname[GADGET_TYPE_GAS]);
          if(ADaPS_exist(plist->data,var_name)){
            calc_mean(  ADaPS_fetch(plist->data,var_name),&mean,   n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);
            calc_median(ADaPS_fetch(plist->data,var_name),&median, n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);
            calc_min(   ADaPS_fetch(plist->data,var_name),&min,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);
            calc_max(   ADaPS_fetch(plist->data,var_name),&max,    n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);
            calc_stddev(ADaPS_fetch(plist->data,var_name),&std_dev,n_of_type[i],SID_REAL,CALC_MODE_RETURN_DOUBLE);
            fprintf(stderr,"%11s  %11.3le  %11.3le  %11.3le  %11.3le  %11.3le K\n",var_name,min,max,mean,median,std_dev);
          }
        }

        /*************/
        /* Write ids */
        /*************/
        for(i=0;i<plist->n_species;i++) {
          if(flag_used[i]){
            sprintf(var_name,"id_%s",pname[i]);
            if(ADaPS_exist(plist->data,var_name)){
              min_i=((size_t *)ADaPS_fetch(plist->data,var_name))[0];
              max_i=((size_t *)ADaPS_fetch(plist->data,var_name))[0];
              for(j=1;j<n_of_type[i];j++){
                min_i=(size_t)MIN(min_i,((size_t *)ADaPS_fetch(plist->data,var_name))[j]);
                max_i=(size_t)MAX(max_i,((size_t *)ADaPS_fetch(plist->data,var_name))[j]);
              }
              fprintf(stderr,"%11s  %11lld  %11lld\n",var_name,min_i,max_i);
            }
          }
        }
      }
      else
        fprintf(stderr,"NO PARTICLES TO ANALYZE!\n");
    }
    fflush(stderr);
    SID_Barrier(SID.COMM_WORLD);
  }
};
