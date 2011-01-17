#include <stdio.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void set_sph_kernel(plist_info *plist,int mode){
  int          i_table;
  int          j_table;
  int          n_temp;
  double      *kernel_radius;
  double      *kernel_radius_temp;
  double      *kernel_table_3d;
  double      *kernel_table_2d;
  double       kernel_table_2d_average;
  interp_info *interp;
  int          flag_compute_2d;
  int          flag_have_3d;
  int          flag_have_2d;
  int          flag_skip;
  
  // Determine if the kernel has already been computed
  flag_compute_2d=check_mode_for_flag(mode,SPH_KERNEL_2D);
  if(!ADaPS_exist(plist->data,"sph_kernel_3d") ||
     (flag_compute_2d && !ADaPS_exist(plist->data,"sph_kernel_2d")))
    flag_skip=FALSE;
  else
    flag_skip=TRUE;

  if(!flag_skip){
    SID_log("Computing SPH kernel...",SID_LOG_OPEN);

    // Initialize kernel arrays
    kernel_radius  =(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
    kernel_table_3d=(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
    for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
      kernel_radius[i_table]  =(double)i_table/(double)N_KERNEL_TABLE;
      kernel_table_3d[i_table]=0.;
    }

    // Set gadget kernel
    if(check_mode_for_flag(mode,SPH_KERNEL_GADGET)){
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
        if(kernel_radius[i_table]<=0.5)
          kernel_table_3d[i_table]= (8./PI)*(1.-6.*kernel_radius[i_table]*kernel_radius[i_table]*(1.-kernel_radius[i_table]));
        else if(kernel_radius[i_table]<=1.)
          kernel_table_3d[i_table]=(16./PI)*(1.-kernel_radius[i_table])*(1.-kernel_radius[i_table])*(1.-kernel_radius[i_table]);
        else
          kernel_table_3d[i_table]=0.;
      }
      ADaPS_store(&(plist->data),(void *)kernel_radius,  "sph_kernel_radius",ADaPS_DEFAULT);
      ADaPS_store(&(plist->data),(void *)kernel_table_3d,"sph_kernel_3d",    ADaPS_DEFAULT);
    }
    else
      SID_trap_error("Unknown kernel type in set_sph_kernel!",ERROR_LOGIC);

    // Integrate to form a projected line-of-sight kernel (if requested)
    if(check_mode_for_flag(mode,SPH_KERNEL_2D)){
      kernel_table_2d   =(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
      kernel_radius_temp=(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
      kernel_table_2d_average=0.;
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
        for(j_table=i_table,n_temp=0;j_table<=N_KERNEL_TABLE;j_table++,n_temp++)
          kernel_radius_temp[n_temp]=sqrt(kernel_radius[j_table]*kernel_radius[j_table]-kernel_radius[i_table]*kernel_radius[i_table]);
        if(n_temp<=1)
          kernel_table_2d[i_table]=0.;
        else if(n_temp<=3){
          init_interpolate(kernel_radius_temp,&(kernel_table_3d[i_table]),(size_t)n_temp,gsl_interp_linear,&interp);
          kernel_table_2d[i_table]=2.*interpolate_integral(interp,0.,kernel_radius_temp[n_temp-1]);
          free_interpolate(&interp);
        }
        else{
          init_interpolate(kernel_radius_temp,&(kernel_table_3d[i_table]),(size_t)n_temp,gsl_interp_cspline,&interp);
          kernel_table_2d[i_table]=2.*interpolate_integral(interp,0.,kernel_radius_temp[n_temp-1]);
          free_interpolate(&interp);
        }
        if(i_table==1)
          kernel_table_2d_average+=PI*kernel_radius[i_table]*kernel_radius[i_table]*0.5*(kernel_table_2d[i_table-1]+kernel_table_2d[i_table]);
        else if(i_table>1)
          kernel_table_2d_average+=TWO_PI*kernel_radius[i_table]*(kernel_radius[i_table]-kernel_radius[i_table-1])*0.5*(kernel_table_2d[i_table-1]+kernel_table_2d[i_table]);
      }
      SID_free((void **)&kernel_radius_temp);
      ADaPS_store(&(plist->data),(void *)kernel_table_2d,           "sph_kernel_2d",   ADaPS_DEFAULT);
      ADaPS_store(&(plist->data),(void *)(&kernel_table_2d_average),"sph_kernel_2d_dA",ADaPS_SCALAR_DOUBLE);
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
}

