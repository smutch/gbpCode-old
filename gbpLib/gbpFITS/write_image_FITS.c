#include <stdio.h>
#include <gbpSID.h>
#include <gbpFITS.h>

int write_image_FITS(void *image,SID_Datatype dtype,int n_D,int *D_in,const char *filename,const char *ext_name){
  fitsfile *fp;
  FILE     *fp_test;
  int   naxis=2;
  int   fpixel=1;
  int   anynull;
  int   i,j;
  long  naxes[2];
  long *D;
  char  keyname[50];
  char  error_msg[80];
  int   status=0;
  long  n_pixels=1;
  int   i_D;
  int   bit_pix;

  // Delete the file if it already exists
  if(fp_test=fopen(filename,"r")) {
    fclose(fp_test);
    remove(filename);
  }

  // Fits files store data in FORTRAN order, so we need to perform a transpose
  transpose_array(image,dtype,n_D,D_in);

  // Convert the array dimensions to a type-long array
  D=(long *)SID_malloc(sizeof(long)*n_D);
  for(i_D=0;i_D<n_D;i_D++)
    D[i_D]=(long)(D_in[i_D]);

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
    bit_pix=FLOAT_IMG;
  }
  // Creeate FLOAT image
  else if(dtype==SID_DOUBLE){
    fits_create_img(fp,DOUBLE_IMG,n_D,D,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS create image: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
    }
    fits_write_img(fp,TDOUBLE,fpixel,n_pixels,image,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS write image: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
    }
    bit_pix=DOUBLE_IMG;
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
    bit_pix=LONG_IMG;
  }

  // Throw error if the given dtype is not supported
  else
    SID_trap_error("Unsupported datatype in write_image_FITS",ERROR_LOGIC);

  // Set extension name
  if(ext_name!=NULL)
    fits_write_key(fp,TSTRING,"EXTNAME",ext_name,NULL,&status);

  // Close the file
  fits_close_file(fp,&status);

  // Fits files store data in FORTRAN order, so we need to perform a transpose
  transpose_array(image,dtype,n_D,D_in);

  SID_free(SID_FARG D);
  return(status);
}

