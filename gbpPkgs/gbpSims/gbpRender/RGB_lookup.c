#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpSPH.h>
#include <gbpCosmo.h>
#include <gbpRender.h>

double RGB_lookup(render_info *render,char colour,int channel){
   return(render->colour_f_RGB[colour][channel]);
}

