#include <math.h>
#include <common.h>

double calc_sep_periodic(double x_1,
		         double y_1,
                         double z_1,
                         double x_2,
                         double y_2,
                         double z_2,
                         double box_size){
  int    i_periodic;
  double x_p,y_p,z_p;
  double r_val;
  r_val=
    (x_2-x_1)*(x_2-x_1)+
    (y_2-y_1)*(y_2-y_1)+
    (z_2-z_1)*(z_2-z_1);

  for(i_periodic=0;i_periodic<26;i_periodic++){
    switch(i_periodic){
    case 0:
      x_p=x_1-box_size;
      y_p=y_1;
      z_p=z_1;
      break;
    case 1:
      x_p=x_1+box_size;
      y_p=y_1;
      z_p=z_1;
      break;
    case 2:
      x_p=x_1;
      y_p=y_1-box_size;
      z_p=z_1;
      break;
    case 3:
      x_p=x_1;
      y_p=y_1+box_size;
      z_p=z_1;
      break;
    case 4:
      x_p=x_1;
      y_p=y_1;
      z_p=z_1-box_size;
      break;
    case 5:
      x_p=x_1;
      y_p=y_1;
      z_p=z_1+box_size;
      break;
    case 6:
      x_p=x_1-box_size;
      y_p=y_1-box_size;
      z_p=z_1;
      break;
    case 7:
      x_p=x_1+box_size;
      y_p=y_1-box_size;
      z_p=z_1;
      break;
    case 8:
      x_p=x_1-box_size;
      y_p=y_1+box_size;
      z_p=z_1;
      break;
    case 9:
      x_p=x_1+box_size;
      y_p=y_1+box_size;
      z_p=z_1;
      break;
    case 10:
      x_p=x_1-box_size;
      y_p=y_1;
      z_p=z_1-box_size;
      break;
    case 11:
      x_p=x_1+box_size;
      y_p=y_1;
      z_p=z_1-box_size;
      break;
    case 12:
      x_p=x_1-box_size;
      y_p=y_1;
      z_p=z_1+box_size;
      break;
    case 13:
      x_p=x_1+box_size;
      y_p=y_1;
      z_p=z_1+box_size;
      break;
    case 14:
      x_p=x_1;
      y_p=y_1-box_size;
      z_p=z_1-box_size;
      break;
    case 15:
      x_p=x_1;
      y_p=y_1+box_size;
      z_p=z_1-box_size;
      break;
    case 16:
      x_p=x_1;
      y_p=y_1-box_size;
      z_p=z_1+box_size;
      break;
    case 17:
      x_p=x_1;
      y_p=y_1+box_size;
      z_p=z_1+box_size;
      break;
    case 18:
      x_p=x_1;
      y_p=y_1-box_size;
      z_p=z_1-box_size;
      break;
    case 19:
      x_p=x_1-box_size;
      y_p=y_1+box_size;
      z_p=z_1-box_size;
      break;
    case 20:
      x_p=x_1-box_size;
      y_p=y_1-box_size;
      z_p=z_1+box_size;
      break;
    case 21:
      x_p=x_1-box_size;
      y_p=y_1+box_size;
      z_p=z_1+box_size;
      break;
    case 22:
      x_p=x_1+box_size;
      y_p=y_1-box_size;
      z_p=z_1-box_size;
      break;
    case 23:
      x_p=x_1+box_size;
      y_p=y_1+box_size;
      z_p=z_1-box_size;
      break;
    case 24:
      x_p=x_1+box_size;
      y_p=y_1-box_size;
      z_p=z_1+box_size;
      break;
    case 25:
      x_p=x_1+box_size;
      y_p=y_1+box_size;
      z_p=z_1+box_size;
      break;
    }
    r_val=MIN(r_val,
	      (x_2-x_p)*(x_2-x_p)+
	      (y_2-y_p)*(y_2-y_p)+
	      (z_2-z_p)*(z_2-z_p));
  }

  return(sqrt(r_val));
}
