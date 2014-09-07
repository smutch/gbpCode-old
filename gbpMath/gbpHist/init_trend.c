#include <stdio.h>
#include <gbpLib.h>
#include <gbpHist.h>

void init_trend(trend_info **trend){
   (*trend)                  =(trend_info *)SID_malloc(sizeof(trend_info));
   (*trend)->ordinate        =NULL;
   (*trend)->coordinate_first=NULL;
   (*trend)->coordinate_last =NULL;
}

