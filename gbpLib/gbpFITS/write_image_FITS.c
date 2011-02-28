#include <gbpSID.h>
#include <gbpFITS.h>

int write_image_FITS(void *image,SID_Datatype dtype,int n_D,long *D,char *filename){
  fitsfile *fp;
  FILE     *fp_test;
  int   naxis=2;
  int   fpixel=1;
  int   anynull;
  int   i,j;
  long  naxes[2];
  char  keyname[50];
  char  error_msg[80];
  int   status=0;
  long  n_pixels=1;
  int   i_D;

  // Delete the file if it already exists
  if(fp_test=fopen(filename,"r")) {
    close(fp);
    remove(filename);
  }

  // Create the file
  fits_create_file(&fp,filename,&status);
  if(status){
    ffgmsg(error_msg);
    SID_trap_error("FITS create file: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
  }

  // Count the number of pixels
  for(i_D=0;i_D<n_D;i_D++)
    n_pixels*=D[i_D];

  // Creeate FLOAT image
  if(dtype==SID_FLOAT){
    fits_create_img(fp,FLOAT_IMG,n_D,D,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS create image: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
    }
    fits_write_img(fp,TFLOAT,fpixel,n_pixels,image,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS write image: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
    }
  }

  // Creeate INTEGER image
  else if(dtype==SID_INT){
    fits_create_img(fp,LONG_IMG,n_D,D,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS create image: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
    }
    fits_write_img(fp,TINT,fpixel,n_pixels,image,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS write image: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
    }
  }

  // Throw error if the given dtype is not supported
  else
    SID_trap_error("Unsupported datatype {%d} in write_image_FITS",ERROR_LOGIC,(int)dtype);


  // Close the file
  fits_close_file(fp,&status);

}

