#include <common.h>
#include <math.h>

// Code from Chris Blake for computing the angle (in degrees) between
//   two points on the sky, given by RA and DEC, expressed in degrees
double calc_angular_separation(double ras1deg,
                               double dec1deg,
                               double ras2deg,
                               double dec2deg){
  double fact,dec1rad,ras1rad,dec2rad,ras2rad;
  double tempreal1,tempreal2,tempreal,distrad;

  fact      = atan(1.0)/45.;
  dec1rad   = dec1deg*fact;
  ras1rad   = ras1deg*fact;
  dec2rad   = dec2deg*fact;
  ras2rad   = ras2deg*fact;
  tempreal1 = sin(dec1rad)*sin(dec2rad);
  tempreal2 = cos(dec1rad)*cos(dec2rad)*cos(ras1rad-ras2rad);
  tempreal  = tempreal1 + tempreal2;
  distrad   = acos(tempreal);

  return(distrad/fact);
}
