#include <gbpLib.h>
#include <gbpFITS.h>

int read_image_FITS(void **image,SID_Datatype *dtype,int *n_D,int **D,char *filename){
  fitsfile *fp;
  int   naxis=2;
  int   fpixel=1;
  int   anynull;
  int   i,j;
  long  naxes[2];
  char  keyname[50];
  char  error_msg[80];
  int   status=0;
  int   n_pixels=0;
  int   i_D;
  int   dtype_fits;
  long  D_in;

  SID_log("Reading FITS image {%s}...",SID_LOG_OPEN,filename);

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
  SID_log("NAXIS =%d",SID_LOG_COMMENT,(*n_D));

  // Read the image dimensions
  (*D)=(int *)SID_malloc(sizeof(int)*(*n_D));
  for(i_D=0;i_D<(*n_D);i_D++){
    sprintf(keyname,"NAXIS%d",i_D+1);
    fits_read_key(fp,TLONG,keyname,&D_in,NULL,&status);
    (*D)[i_D]=D_in;
    if(status){
      ffgmsg(error_msg);
      SID_trap_error("FITS read {%s}: filename={%s} status=%d message={%s}",ERROR_IO_OPEN,keyname,filename,status,error_msg);
    }
    if(i_D==0)
      n_pixels=(*D)[i_D];
    else
      n_pixels*=(*D)[i_D];
    SID_log("NAXIS%d=%d",SID_LOG_COMMENT,i_D,(*D)[i_D]);
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
      SID_log("dtype =INTEGER",SID_LOG_COMMENT);
      break;
    case 32:
      (*dtype)=SID_LONG_LONG;
      (*image)=SID_malloc(sizeof(long long)*n_pixels);
      fits_read_img(fp,TLONG,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      SID_log("dtype =LONG LONG",SID_LOG_COMMENT);
      break;
    case -32:
      (*dtype)=SID_FLOAT;
      (*image)=SID_malloc(sizeof(float)*n_pixels);
      fits_read_img(fp,TFLOAT,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      SID_log("dtype =FLOAT",SID_LOG_COMMENT);
      break;
    case -64:
      (*dtype)=SID_DOUBLE;
      (*image)=SID_malloc(sizeof(double)*n_pixels);
      fits_read_img(fp,TDOUBLE,fpixel,n_pixels,0,(*image),&anynull,&status);
      if(status){
        ffgmsg(error_msg);
        SID_trap_error("FITS image: filename={%s} dtype_fits=%d status=%d message={%s}",ERROR_IO_OPEN,filename,dtype_fits,status,error_msg);
      }
      SID_log("dtype =DOUBLE",SID_LOG_COMMENT);
      break;
    default:
      SID_trap_error("Unsupported datatype {%d} in write_image_FITS",ERROR_LOGIC,dtype_fits);
      break;
  }

  // Close the file
  fits_close_file(fp,&status);
  SID_log("Done.",SID_LOG_CLOSE);

}

