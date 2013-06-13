#include <gbpLib.h>
#include <gbpMisc.h>
#include <gbpSort.h>
#include <string.h>
#include <math.h>

double calc_weight_local(double *W,int i);
double calc_weight_local(double *W,int i){
   if(W!=NULL)
      return(W[i]);
   else
      return(1.);
}

int compute_centroid3D(double  *W,
                       double  *x_in,
                       double  *y_in,
	               double  *z_in,
                       int      n,
                       double   R_min,
                       double   step,
                       int      convergence_N,
                       int      mode,
                       double  *xcen_out,
	   	       double  *ycen_out,
	   	       double  *zcen_out){
  int n_iterations=0;

  // Parse the mode flags
  int flag_centroid_step;
  int flag_centroid_inplace;
  int flag_return_indices;
  if(check_mode_for_flag(mode,CENTROID3D_MODE_STEP)){
     flag_centroid_step=TRUE;
     if(step<=0.)
        SID_trap_error("Invalid step size selected (%le) for mode %d in compute_centroid3D().",ERROR_LOGIC,step,mode);
  }
  else if(check_mode_for_flag(mode,CENTROID3D_MODE_FACTOR)){
     flag_centroid_step=FALSE;
     if(step>=1. || step<=0.)
        SID_trap_error("Invalid step size selected (%le) for mode %d in compute_centroid3D().",ERROR_LOGIC,step,mode);
  }
  else
     SID_trap_error("Invalid mode (%d) specified in compute_3dcentroid().",ERROR_LOGIC,mode);
  flag_centroid_inplace=check_mode_for_flag(mode,CENTROID3D_MODE_INPLACE);
  flag_return_indices  =check_mode_for_flag(mode,CENTROID3D_MODE_RETURN_INDICES);

  // Perform centroiding if n>0
  double  xcen;
  double  ycen;
  double  zcen;
  double *x;
  double *y;
  double *z;
  double *R;
  size_t *R_index;
  size_t  i,j;
  xcen=0.0;
  ycen=0.0;
  zcen=0.0;
  if(n>0){
    if(flag_centroid_inplace){
       x=x_in;
       y=y_in;
       z=z_in;
    }
    else{
       x=(double *)SID_malloc(sizeof(double)*n);
       y=(double *)SID_malloc(sizeof(double)*n);
       z=(double *)SID_malloc(sizeof(double)*n);
       memcpy(x,x_in,sizeof(double)*n);
       memcpy(y,y_in,sizeof(double)*n);
       memcpy(z,z_in,sizeof(double)*n);
    }
    R=(double *)SID_malloc(sizeof(double)*n);

    // Compute starting point
    for(i=0;i<n;i++){
       xcen+=x[i];
       ycen+=y[i];
       zcen+=z[i];
    }
    xcen/=(double)n;
    ycen/=(double)n;
    zcen/=(double)n;
    for(i=0;i<n;i++){
       x[i]-=xcen;
       y[i]-=ycen;
       z[i]-=zcen;
    }

    // Compute initial radii
    double R_max;
    int    N;
    for(i=0,N=0;i<n;i++){
      R[i]=pow(x[i]*x[i]+
               y[i]*y[i]+
               z[i]*z[i],0.5);
      N++;
    }

    // Sort initial radii
    merge_sort(R,(size_t)n,&R_index,SID_DOUBLE,SORT_COMPUTE_INDEX,FALSE);
    R_max=R[R_index[n-1]];

    // Initialize convergence criteria and some working variables
    double R_ap;
    double W_tot;
    double x_next;
    double y_next;
    double z_next;
    x_next=0.;
    y_next=0.;
    z_next=0.;
    R_ap  =R_max;

    // Loop until convergence
    int continue_flag;
    if(N>convergence_N && R_ap>R_min)
       continue_flag=TRUE;
    else
       continue_flag=FALSE;
    while(continue_flag){

      // Recentre particles and recompute radii
      if(n_iterations>0){
         for(i=0;i<n;i++){
           x[i]-=x_next;
           y[i]-=y_next;
           z[i]-=z_next;
           R[i] =pow(x[i]*x[i]+
                     y[i]*y[i]+
                     z[i]*z[i],0.5);
         }

         // Sort radii
         SID_free(SID_FARG R_index);
         merge_sort(R,(size_t)n,&R_index,SID_DOUBLE,SORT_COMPUTE_INDEX,FALSE);
      }

      // Test for a new centre
      double weight_i;
      x_next=0.0;
      y_next=0.0;
      z_next=0.0;
      W_tot =0.0;
      N     =0;
      j     =0;
      i     =R_index[j];
      for(;j<n && R[i]<=R_ap ;j++){
        i        =R_index[j];
        weight_i =calc_weight_local(W,i);
        x_next  +=weight_i*x[i];
        y_next  +=weight_i*y[i];
        z_next  +=weight_i*z[i];
        W_tot   +=weight_i;
        N++;
      }
      if(N>=convergence_N && R_ap>R_min){
        x_next/=W_tot;
        y_next/=W_tot;
        z_next/=W_tot;
        xcen  +=x_next;
        ycen  +=y_next;
        zcen  +=z_next;
        if(flag_centroid_step)
          R_ap-=step;
        else
          R_ap*=step;
      }
      else
        continue_flag=FALSE;
      n_iterations++;
    }

    // Clean-up
    if(!flag_centroid_inplace){
       SID_free(SID_FARG W);
       SID_free(SID_FARG x);
       SID_free(SID_FARG y);
       SID_free(SID_FARG z);
    }
    SID_free(SID_FARG R);
    SID_free(SID_FARG R_index);
  }

  // Finalize return values
  (*xcen_out)=xcen;
  (*ycen_out)=ycen;
  (*zcen_out)=zcen;

  return(n_iterations);
}

