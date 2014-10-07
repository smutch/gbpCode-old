#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpParse_core.h>

void float_to_text(float number,int n_sig_fig,char *result){
   int   number0=(int)number;
   float number1=number-(float)number0;
   float number2=number1;
   for(int i_sig_fig=0;i_sig_fig<n_sig_fig;i_sig_fig++)
      number2*=10;
   int number3=(int)number2;
   sprintf(result,"%dpt%d",number0,number3);
}

