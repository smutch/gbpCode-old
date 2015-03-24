#include <stdio.h>
#include <gbpLib.h>
#include <gbpRender.h>

void set_sph_kernel(double **kernel_radius,
                    double **kernel_table_3d,
                    double **kernel_table_2d,
                    double  *kernel_table_2d_average,
                    int      mode){
  int          i_table;
  int          j_table;
  int          n_temp;
  interp_info *interp;
  int          flag_compute_2d;
  int          flag_have_3d;
  int          flag_have_2d;
  int          flag_skip;
  double      *kernel_radius_temp;
  
  // Determine if the kernel has already been computed
  flag_compute_2d=check_mode_for_flag(mode,SPH_KERNEL_2D);
  if((*kernel_table_3d)==NULL || (flag_compute_2d && (*kernel_table_2d)==NULL))
    flag_skip=FALSE;
  else
    flag_skip=TRUE;

  if(!flag_skip){
    SID_log("Computing SPH kernel...",SID_LOG_OPEN|SID_LOG_TIMER);

    // Initialize kernel arrays
    (*kernel_radius)  =(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
    (*kernel_table_3d)=(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));

    // Set gadget kernel
    if(check_mode_for_flag(mode,SPH_KERNEL_GADGET)){
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
        (*kernel_radius)[i_table]  =(double)i_table/(double)N_KERNEL_TABLE;
        (*kernel_table_3d)[i_table]=0.;
      }
      // From subfind ...
      //#define  NUMDIMS 3      //!< For 3D-normalized kernel 
      //#define  KERNEL_COEFF_1  2.546479089470 //!< Coefficients for SPH spline kernel and its derivative 
      //#define  KERNEL_COEFF_2  15.278874536822
      //#define  KERNEL_COEFF_3  45.836623610466
      //#define  KERNEL_COEFF_4  30.557749073644
      //#define  KERNEL_COEFF_5  5.092958178941
      //#define  KERNEL_COEFF_6  (-15.278874536822)
      //#define  NORM_COEFF      4.188790204786 //!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 
      //u = r * hinv;
      //if(u < 0.5)
      //  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
      //else
      //  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
      //mass_j = P_Mass;
      //rho += (mass_j * wk);
      // ... but I code things according to the paper ...
      // n.b.: m_p*(1/h)^3 must be applied later
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
        if((*kernel_radius)[i_table]<=0.5)
          (*kernel_table_3d)[i_table]= (8./PI)*(1.-6.*(*kernel_radius)[i_table]*(*kernel_radius)[i_table]*(1.-(*kernel_radius)[i_table]));
        else if((*kernel_radius)[i_table]<=1.)
          (*kernel_table_3d)[i_table]=(16./PI)*(1.-(*kernel_radius)[i_table])*(1.-(*kernel_radius)[i_table])*(1.-(*kernel_radius)[i_table]);
        else
          (*kernel_table_3d)[i_table]=0.;
      }
    }
    else if(check_mode_for_flag(mode,SPH_KERNEL_GAUSSIAN)){
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
        (*kernel_radius)[i_table]  =3.*(double)i_table/(double)N_KERNEL_TABLE;
        (*kernel_table_3d)[i_table]=0.;
      }
      // Radius is in units of standard deviation
      double norm=1./sqrt(TWO_PI);
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++)
         (*kernel_table_3d)[i_table]=norm*exp(-0.5*(*kernel_radius)[i_table]*(*kernel_radius)[i_table]);
    }
    else
      SID_trap_error("Unknown kernel type in set_sph_kernel!",ERROR_LOGIC);

    // Integrate to form a projected line-of-sight kernel (if requested)
    if(check_mode_for_flag(mode,SPH_KERNEL_2D)){
      (*kernel_table_2d)=(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
      kernel_radius_temp=(double *)SID_malloc(sizeof(double)*(N_KERNEL_TABLE+1));
      (*kernel_table_2d_average)=0.;
      for(i_table=0;i_table<=N_KERNEL_TABLE;i_table++){
        for(j_table=i_table,n_temp=0;j_table<=N_KERNEL_TABLE;j_table++,n_temp++)
          kernel_radius_temp[n_temp]=sqrt((*kernel_radius)[j_table]*(*kernel_radius)[j_table]-(*kernel_radius)[i_table]*(*kernel_radius)[i_table]);
        if(n_temp<=1)
          (*kernel_table_2d)[i_table]=0.;
        else if(n_temp<=3){
          init_interpolate(kernel_radius_temp,&((*kernel_table_3d)[i_table]),(size_t)n_temp,gsl_interp_linear,&interp);
          (*kernel_table_2d)[i_table]=2.*interpolate_integral(interp,0.,kernel_radius_temp[n_temp-1]);
          free_interpolate(SID_FARG interp,NULL);
        }
        else{
          init_interpolate(kernel_radius_temp,&((*kernel_table_3d)[i_table]),(size_t)n_temp,gsl_interp_cspline,&interp);
          (*kernel_table_2d)[i_table]=2.*interpolate_integral(interp,0.,kernel_radius_temp[n_temp-1]);
          free_interpolate(SID_FARG interp,NULL);
        }
        if(i_table==1)
          (*kernel_table_2d_average)+=PI*(*kernel_radius)[i_table]*(*kernel_radius)[i_table]*0.5*((*kernel_table_2d)[i_table-1]+(*kernel_table_2d)[i_table]);
        else if(i_table>1)
          (*kernel_table_2d_average)+=TWO_PI*(*kernel_radius)[i_table]*((*kernel_radius)[i_table]-(*kernel_radius)[i_table-1])*0.5*((*kernel_table_2d)[i_table-1]+(*kernel_table_2d)[i_table]);
      }
      SID_free(SID_FARG kernel_radius_temp);
    }

    SID_log("Done.",SID_LOG_CLOSE);
  }
}

