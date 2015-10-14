#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  int         i;
  plist_info  plist;
  int         snapshot;
  char        filename_in[256];
  char        filename_out[256];
  int         flag_mode;
  int         n_input_vals;
  size_t      n_marked;
  double     *input_vals;
  double      h_Hubble;

  SID_init(&argc,&argv,NULL,NULL);

  // Parse command line
  if(argc<=1){
    SID_log_error("\nsyntax: %s filename_in_root snapshot filename_out [box/sphere] [x_min,xmax,y_min,y_max,z_min,z_max/x_cen,y_cen,z_cen,radius]",ERROR_SYNTAX,argv[0]);
    SID_trap_error("------\n",ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in, argv[1]);
    snapshot=atoi(argv[2]);
    strcpy(filename_out,argv[3]);
    if(!strcmp(argv[4],"box"))
      flag_mode=VOLUME_BOX;
    else if(!strcmp(argv[4],"sphere"))
      flag_mode=VOLUME_SPHERE;
    else
      SID_trap_error("Valid modes are \"box\" or \"sphere\".",ERROR_SYNTAX);
    if(flag_mode==VOLUME_BOX)
      n_input_vals=6;
    else if(flag_mode==VOLUME_SPHERE)
      n_input_vals=4;
  }

  SID_log("Marking particles in GADGET file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);

  // Initialize input parameters
  input_vals=(double *)SID_malloc(sizeof(double)*n_input_vals);
  for(i=0;i<n_input_vals;i++)
    input_vals[i]=(double)atof(argv[5+i]);

  // Report selection volume
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double z_min;
  double z_max;
  switch(flag_mode){
    case VOLUME_BOX:
      SID_log("Selection Box:",SID_LOG_OPEN);
      SID_log("x_min: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[0]);
      SID_log("x_max: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[1]);
      SID_log("y_min: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[2]);
      SID_log("y_max: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[3]);
      SID_log("z_min: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[4]);
      SID_log("z_max: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[5]);
      x_min=input_vals[0];
      x_max=input_vals[1];
      y_min=input_vals[2];
      y_max=input_vals[3];
      z_min=input_vals[4];
      z_max=input_vals[5];
      break;
    case VOLUME_SPHERE:
      SID_log("Selection Sphere:",SID_LOG_OPEN);
      SID_log("x_centre: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[0]);
      SID_log("y_centre: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[1]);
      SID_log("z_centre: %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[2]);
      SID_log("radius:   %10.3lf [Mpc/h]",SID_LOG_COMMENT,input_vals[3]);
      x_min=input_vals[0]-input_vals[3];
      x_max=input_vals[0]+input_vals[3];
      y_min=input_vals[1]-input_vals[3];
      y_max=input_vals[1]+input_vals[3];
      z_min=input_vals[2]-input_vals[3];
      z_max=input_vals[2]+input_vals[3];
      break;
    default:
      SID_trap_error("Unknown selection mode",ERROR_LOGIC);
      break;
  }
  SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);

  // Read GADGET file
  SID_log("Reading GADGET file...",SID_LOG_OPEN|SID_LOG_TIMER);
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_gadget_binary_header(filename_in,snapshot,&plist);
  h_Hubble=((double *)ADaPS_fetch(plist.data,"h_Hubble"))[0];
  x_min*=M_PER_MPC/h_Hubble;
  x_max*=M_PER_MPC/h_Hubble;
  y_min*=M_PER_MPC/h_Hubble;
  y_max*=M_PER_MPC/h_Hubble;
  z_min*=M_PER_MPC/h_Hubble;
  z_max*=M_PER_MPC/h_Hubble;
  ADaPS_store(&(plist.data),&x_min,"rank_x_min",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),&x_max,"rank_x_max",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),&y_min,"rank_y_min",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),&y_max,"rank_y_max",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),&z_min,"rank_z_min",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),&z_max,"rank_z_max",ADaPS_SCALAR_DOUBLE);
  read_gadget_binary(filename_in,snapshot,&plist,READ_GADGET_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);

  // Convert units
  for(i=0;i<n_input_vals;i++)
    input_vals[i]*=M_PER_MPC/h_Hubble;

  // Mark particles
  SID_log("Marking particles...",SID_LOG_OPEN|SID_LOG_TIMER);
  n_marked=mark_particles(&plist,flag_mode|MARK_INIT,input_vals,"mark");
  SID_log("%lld particles marked...Done.",SID_LOG_CLOSE,n_marked);

  // Write mark file
  SID_log("Writing mark file...",SID_LOG_OPEN|SID_LOG_TIMER);
  write_mark_file(&plist,"mark",filename_out);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_plist(&plist);
  if(n_input_vals>0)
    SID_free((void **)&input_vals);
  SID_log("Done.",SID_LOG_CLOSE);

  return(ERROR_NONE);
}
