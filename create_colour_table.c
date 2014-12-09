#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpRender.h>

int integer_gaussian_local(int i_colour,float f_amplitude,float f_centre,float f_width,int rgb_max,int n_colours);
int integer_gaussian_local(int i_colour,float f_amplitude,float f_centre,float f_width,int rgb_max,int n_colours){
   return((int)(f_amplitude*rgb_max*exp(-pow(((float)(i_colour)-f_centre*(float)n_colours)/(float)(f_width*(float)n_colours),2.))));
}

/************************************************************/
/* If colourmapselect<0 then invert                         */
/************************************************************/
void create_colour_table(int     colourmapselect,
                         int     n_colours,
                         int  ***rgb){
  int   i,j;
  int  *rgb_temp;
  int   colourmap;
  int   flag_invert;
  int   rgb_max=255;

  // Interpret selection
  if(colourmapselect<0){
    colourmap  =-1*colourmapselect;
    flag_invert=TRUE;
  }
  else{
    colourmap  =colourmapselect;
    flag_invert=FALSE;
  }

  // Initialize array
  (*rgb)   =(int **)SID_malloc(sizeof(int *)*3);
  (*rgb)[0]=(int  *)SID_malloc(sizeof(int)*n_colours);
  (*rgb)[1]=(int  *)SID_malloc(sizeof(int)*n_colours);
  (*rgb)[2]=(int  *)SID_malloc(sizeof(int)*n_colours);
  for(i=0;i<n_colours;i++){
    (*rgb)[0][i]=0;
    (*rgb)[1][i]=0;
    (*rgb)[2][i]=0;
  }

  // Create colour table
  switch(colourmap){
  case 1:
    for(i=0;i<n_colours;i++){
      (*rgb)[0][i]=(int)((float)rgb_max*(float)i/(float)(n_colours-1));
      (*rgb)[1][i]=(*rgb)[0][i];
      (*rgb)[2][i]=(*rgb)[0][i];
    }
    break;
  case 2:
    for(i=0;i<n_colours;i++){
      (*rgb)[0][i]=(int)((float)rgb_max*(float)i/(float)(n_colours-1));
      (*rgb)[1][i]=0;
      (*rgb)[2][i]=0;
    }
    break;
  case 3:
    for(i=0;i<n_colours;i++){
      (*rgb)[1][i]=(int)((float)rgb_max*(float)i/(float)(n_colours-1));
      (*rgb)[0][i]=0;
      (*rgb)[2][i]=0;
    }
    break;
  case 4:
    for(i=0;i<n_colours;i++){
      (*rgb)[2][i]=(int)((float)rgb_max*(float)i/(float)(n_colours-1));
      (*rgb)[0][i]=0;
      (*rgb)[1][i]=0;
    }
    break;
  case 5:
    i=0;
    for(j=0;j<n_colours;j+=3)
      (*rgb)[2][i++]=j;
    for(j=n_colours-1;j>0 && i<n_colours;j-=3) {
      (*rgb)[0][i]=(*rgb)[0][i-1]+3;
      (*rgb)[2][i]=j;
      i++;
    }
    while(i<n_colours) {
      (*rgb)[1][i]=(*rgb)[1][i-1]+3;
      (*rgb)[0][i]=(*rgb)[0][i-1]+3;
      i++;
    }
    break;
  case 6:
    i=1;
    for(j=1;j<n_colours && i<n_colours;j++)
      (*rgb)[2][i++]=rgb_max-(int)(1.8*(float)j);
    i=1;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=1;j<n_colours && i<n_colours;j++)
      (*rgb)[0][i++]=(int)((float)j*2.0);
    i=(int)((double)(n_colours)/2);
    j=i;
    while(i<n_colours) {
      (*rgb)[1][i]=(int)((float)(i-j)*1.7);
      i++;
    }
    break;
  case 7:
    for(i=0;i<n_colours;i++){      
      (*rgb)[0][i]=(int)(rgb_max*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/2.25),2.)));
      if(i<n_colours/4)
        (*rgb)[1][i]=(int)(rgb_max*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/2),2.)));
      else
        (*rgb)[1][i]=rgb_max;
      if(i<n_colours/4)
        (*rgb)[2][i]=MAX((int)(200.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/2),2.))),
                         (int)(200.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/7),2.))));
      else
        (*rgb)[2][i]=MAX((int)(200.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/4),2.))),
                         (int)(200.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/7),2.))));
    }
    break;
  case 8:
    for(i=0;i<n_colours;i++){
      if(i<3*n_colours/4)
        (*rgb)[0][i]=(int)((double)rgb_max*exp(-pow((float)(i-3*(int)rgb_max/4)/(float)(rgb_max/5),2.)));
      else
        (*rgb)[0][i]=rgb_max;
      if(i<n_colours/3)
        (*rgb)[1][i]=(int)((double)rgb_max*exp(-pow(((float)(i-(int)rgb_max/3))/(float)(rgb_max/4),2.)));
      else
        (*rgb)[1][i]=rgb_max;
      if(i<n_colours/4)
        (*rgb)[2][i]=rgb_max;
      else
        (*rgb)[2][i]=MAX((int)((double)rgb_max*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/5),2.))),
                         (int)((double)rgb_max*exp(-pow((float)(i-(int)rgb_max)  /(float)(rgb_max/8),2.))));
    }
    break;
  case 9:
    for(i=0;i<n_colours;i++){
      (*rgb)[0][i]=(int)(rgb_max*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/2.25),2.)));
      if(i<n_colours/4)
        (*rgb)[1][i]=(int)(rgb_max*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/4),2.)));
      else
        (*rgb)[1][i]=rgb_max;
      if(i<n_colours/4)
        (*rgb)[2][i]=MAX((int)(200.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/4),2.))),
                         (int)(200.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/7),2.))));
      else
        (*rgb)[2][i]=MAX((int)(200.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/4),2.))),
                         (int)(200.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/7),2.))));
    }
    break;
  case 10:
    i=1;
    for(j=0;j<n_colours && i<n_colours;j++)
      (*rgb)[2][i++]=rgb_max;
    i=1;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=0;j<n_colours && i<n_colours;j++)
      (*rgb)[0][i++]=(int)((float)j*2.0);
    i=(int)((double)(n_colours)/2);
    j=i;
    while(i<n_colours) 
      (*rgb)[1][i++]=(int)((float)(i-j)*1.7);
    break;
  case 11:
    i=1;
    for(j=0;j<n_colours && i<n_colours;j++)
      (*rgb)[2][i++]=rgb_max-(int)(1.5*(float)j);
    i=1;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=0;j<n_colours && i<n_colours;j++)
      (*rgb)[0][i++]=rgb_max;
    i=(int)((double)(rgb_max)/2.5);
    j=i;
    while(i<n_colours) 
      (*rgb)[1][i++]=(int)((float)(i-j)*1.7);
    break;
  case 12:
    i=0;
    for(j=0;j<n_colours && i<n_colours;j++)
      (*rgb)[2][i++]=rgb_max-(int)(1.5*(float)j);
    i=0;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=0;j<n_colours && i<n_colours;j++)
      (*rgb)[0][i++]=rgb_max;
    i=(int)((double)(rgb_max)/2.5);
    j=i;
    while(i<n_colours)
      (*rgb)[1][i++]=(int)((float)(i-j)*1.7);
    break;
  //integer_gaussian_local(int i_colour,float f_amplitude,float f_centre,float f_width,int rgb_max,int n_colours);
  case 13:
    for(i=0;i<n_colours;i++){
      //if(i<(n_colours/10))
      //   (*rgb)[0][i]=0.8*rgb_max;
      //else
      //   (*rgb)[0][i]=integer_gaussian_local(i,0.8,0.1,0.50,rgb_max,n_colours);
      //(*rgb)[1][i]   =integer_gaussian_local(i,0.8,0.0,0.25,rgb_max,n_colours);
      ////(*rgb)[2][i]   =integer_gaussian_local(i,1.0,0.0,0.15,rgb_max,n_colours);
      //(*rgb)[0][i]  +=integer_gaussian_local(i,0.1,1.0,0.20,rgb_max,n_colours);
      //(*rgb)[1][i]  +=integer_gaussian_local(i,0.3,0.6,0.50,rgb_max,n_colours);
      //(*rgb)[2][i]  +=integer_gaussian_local(i,0.5,1.0,0.20,rgb_max,n_colours);

      if(i<(n_colours/10))
         (*rgb)[0][i]=1.*rgb_max;
      else
         (*rgb)[0][i]=integer_gaussian_local(i,1.0,0.1,0.50,rgb_max,n_colours);
      (*rgb)[1][i]   =integer_gaussian_local(i,1.0,0.0,0.25,rgb_max,n_colours);
      //(*rgb)[2][i]   =integer_gaussian_local(i,1.0,0.0,0.15,rgb_max,n_colours);
      (*rgb)[0][i]  +=integer_gaussian_local(i,0.4,1.0,0.20,rgb_max,n_colours);
      (*rgb)[1][i]  +=integer_gaussian_local(i,0.6,0.6,0.50,rgb_max,n_colours);
      (*rgb)[2][i]  +=integer_gaussian_local(i,1.0,1.0,0.20,rgb_max,n_colours);

    }
    break;
  }

  // Make sure everything is within bounds
  for(i=0;i<n_colours;i++) {
    (*rgb)[0][i]=MAX(0,MIN((*rgb)[0][i],n_colours-1));
    (*rgb)[1][i]=MAX(0,MIN((*rgb)[1][i],n_colours-1));
    (*rgb)[2][i]=MAX(0,MIN((*rgb)[2][i],n_colours-1));
  }

  // Invert colourtable (if required)
  if(colourmapselect<0){
    rgb_temp=(int *)SID_malloc(sizeof(int)*(n_colours));
    for(j=0;j<n_colours;j++)
      rgb_temp[j]=(*rgb)[0][n_colours-j-1];
    for(j=0;j<n_colours;j++)
      (*rgb)[0][j]=rgb_temp[j];
    for(j=0;j<n_colours;j++)
      rgb_temp[j]=(*rgb)[1][n_colours-j-1];
    for(j=0;j<n_colours;j++)
      (*rgb)[1][j]=rgb_temp[j];
    for(j=0;j<n_colours;j++)
      rgb_temp[j]=(*rgb)[2][n_colours-j-1];
    for(j=0;j<n_colours;j++)
      (*rgb)[2][j]=rgb_temp[j];
    SID_free(SID_FARG rgb_temp);
  }

}
