#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>
#include <gbpRender.h>

/*********************************************************/
/* |colourmapselect| = 0 -> discrete                     */
/*                     1 -> greyscale                    */
/*                     2 -> blue,red,orange,yellow       */
/*                     3 -> black,blue,red,orange,yellow */
/*                     4 -> cyan->green->yellow->white   */
/*                     5 -> blue->green->yellow->white   */
/*                     6 -> blue,pink,white              */
/*                     7 -> pink,red,yellow              */
/* If colourmapselect<0 then invert                      */
/*********************************************************/
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
  case 3:
    i=1;
    for(j=1;j<n_colours;j++)
      (*rgb)[2][i++]=rgb_max-(int)(1.8*(float)j);
    i=1;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=1;j<n_colours;j++)
      (*rgb)[0][i++]=(int)((float)j*2.0);
    i=(int)((double)(n_colours)/2);
    j=i;
    while(i<n_colours) {
      (*rgb)[1][i]=(int)((float)(i-j)*1.7);
      i++;
    }
    break;
  case 4:
    for(i=0;i<n_colours;i++){      
      (*rgb)[0][i]=(int)(255.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/2.25),2.)));
      if(i<n_colours/4)
        (*rgb)[1][i]=(int)(255.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/2),2.)));
      else
        (*rgb)[1][i]=255;
      if(i<n_colours/4)
        (*rgb)[2][i]=MAX((int)(200.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/2),2.))),
                         (int)(200.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/7),2.))));
      else
        (*rgb)[2][i]=MAX((int)(200.*exp(-pow((float)(i-(int)rgb_max/4)/(float)(rgb_max/4),2.))),
                         (int)(200.*exp(-pow((float)(i-(int)rgb_max)/(float)(rgb_max/7),2.))));
    }
    break;
  case 5:
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
  case 6:
    i=1;
    for(j=0;j<n_colours;j++)
      (*rgb)[2][i++]=rgb_max;
    i=1;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=0;j<n_colours;j++)
      (*rgb)[0][i++]=(int)((float)j*2.0);
    i=(int)((double)(n_colours)/2);
    j=i;
    while(i<n_colours) 
      (*rgb)[1][i++]=(int)((float)(i-j)*1.7);
    break;
  case 7:
    i=1;
    for(j=0;j<n_colours;j++)
      (*rgb)[2][i++]=rgb_max-(int)(1.5*(float)j);
    i=1;
    (*rgb)[0][i]=(int)((double)(rgb_max)/5.0);
    for(j=0;j<n_colours;j++)
      (*rgb)[0][i++]=rgb_max;
    i=(int)((double)(rgb_max)/2.5);
    j=i;
    while(i<n_colours) 
      (*rgb)[1][i++]=(int)((float)(i-j)*1.7);
    break;
  case 2:
  default:
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
