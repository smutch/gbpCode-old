#include <gbpLib.h>
#include <gbpFITS.h>

int read_image_FITS(void **image,SID_Datatype *dtype,int *n_D,long **D,char *filename){
  fitsfile *fp;
  int   naxis=2;
  int   fpixel=1;
  int   anynull;
  int   i,j;
  long  naxes[2];
  char  keyname[50];
  char  error_msg[80];
  int   status;
  int   n_pixels=1;
  int   i_D;
  int   dtype_fits;

  // Create the file
  fits_open_file(&fp,filename,READONLY,&status);	      
  if(status){
    ffgmsg(error_msg);
    SID_trap_error("FITS open file: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,filename,status,error_msg);
  }

  sprintf(keyname,"NAXIS");
  fits_read_key(fp,TINT,keyname,n_D,NULL,&status);
  if(status){
    ffgmsg(error_msg);
    SID_trap_error("FITS read {%s}: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,keyname,filename,status,error_msg);
  }

  // Read the image dimensions
  (*D)=(long *)SID_malloc(sizeof(long)*(*n_D));
  for(i_D=0;i_D<(*n_D);i_D++){
    sprintf(keyname,"NAXIS%d",i_D);
    fits_read_key(fp,TLONG,keyname,&((*D)[i_D]),NULL,&status);
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS read {%s}: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,keyname,filename,status,error_msg);
    }
    n_pixels*=(*D)[i_D];
  }

  // Read the data type
  sprintf(keyname,"BITPIX");
  fits_read_key(fp,TINT,keyname,&dtype_fits,NULL,&status);
  if(status){
    ffgmsg(error_msg);
    SID_trap_error("FITS read {%s}: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,keyname,filename,status,error_msg);
  }
  switch(dtype_fits){
    case 16:
      (*dtype)=SID_INT;
      (*image)=SID_malloc(sizeof(int)*n_pixels);
      fits_read_img(fp,TINT,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      break;
    case 32:
      (*dtype)=SID_LONG_LONG;
      (*image)=SID_malloc(sizeof(long long)*n_pixels);
      fits_read_img(fp,TLONG,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      break;
    case -32:
      (*dtype)=SID_FLOAT;
      (*image)=SID_malloc(sizeof(float)*n_pixels);
      fits_read_img(fp,TFLOAT,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      break;
    case -64:
      (*dtype)=SID_DOUBLE;
      (*image)=SID_malloc(sizeof(double)*n_pixels);
      fits_read_img(fp,TDOUBLE,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      break;
    default:
      SID_trap_error("Unsupported datatype {%d} in write_image_FITS",ERROR_LOGIC,dtype_fits);
      break;
  }

  // Close the file
  fits_close_file(fp,&status);

}

