// Free memory allocated by initialize_field()
#include <gbpLib.h>
#include <gbpPHKs.h>

PHK_t compute_PHK_from_Cartesian(int PHK_bit_size,int n_D,...){
  int      i_d;
  double   x=0.;
  double   y=0.;
  double   z=0.;
  va_list  vargs;
  va_start(vargs,n_D);

  if(n_D>3 || n_D<0)
    SID_trap_error("Invalid number of dimensions (%d) in compute_PHK_from_Cartesian",ERROR_LOGIC,n_D);
  for(i_d=0;i_d<n_D;i_d++){
    switch(i_d){
      case 0:
        x=(double)va_arg(vargs,double);
        if(x==1.)
          x=0.;
        else if(x<0. || x>=1.)
          SID_trap_error("x (%le) is not in the range [0,1) in compute_PHK_from_Cartesian",ERROR_LOGIC,x);
        break;
      case 1:
        y=(double)va_arg(vargs,double);
        if(y==1.)
          y=0.;
        else if(y<0. || y>=1.)
          SID_trap_error("y (%le) is not in the range [0,1) in compute_PHK_from_Cartesian",ERROR_LOGIC,y);
        break;
      case 2:
        z=(double)va_arg(vargs,double);
        if(z==1.)
          z=0.;
        else if(z<0. || z>=1.)
          SID_trap_error("z (%le) is not in the range [0,1) in compute_PHK_from_Cartesian",ERROR_LOGIC,z);
        break;
    }
  }
  va_end(vargs);
  return((PHK_t)sfc_curve_calcKey(SFC_CURVE_HILBERT,x,y,z,PHK_bit_size));
}

