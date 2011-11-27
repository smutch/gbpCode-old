#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
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
  double      h_Hubble=1.; // Don't apply h
  double      x_min,x_max;
  double      y_min,y_max;
  double      z_min,z_max;

  SID_init(&argc,&argv,NULL);

  // Parse command line
  if(argc<=1){
    SID_log_error("\nsyntax: %s filename_in snapshot filename_out [box/sphere] [x_min,xmax,y_min,y_max,z_min,z_max/x_cen,y_cen,z_cen,radius]",ERROR_SYNTAX,argv[0]);
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
    else if(flag_mode=VOLUME_SPHERE)
      n_input_vals=4;
  }

  // Initialize input parameters
  input_vals=(double *)SID_malloc(sizeof(double)*n_input_vals);
  for(i=0;i<n_input_vals;i++)
    input_vals[i]=(double)atof(argv[4+i])*M_PER_MPC/h_Hubble;
  if(flag_mode==VOLUME_BOX){
    x_min=input_vals[0]*M_PER_MPC;
    x_max=input_vals[1]*M_PER_MPC;
    y_min=input_vals[2]*M_PER_MPC;
    y_max=input_vals[3]*M_PER_MPC;
    z_min=input_vals[4]*M_PER_MPC;
    z_max=input_vals[5]*M_PER_MPC;
  }
  else if(flag_mode==VOLUME_SPHERE){
    x_min=(input_vals[0]-input_vals[3])*M_PER_MPC;
    x_max=(input_vals[0]+input_vals[3])*M_PER_MPC;
    y_min=(input_vals[1]-input_vals[3])*M_PER_MPC;
    y_max=(input_vals[1]+input_vals[3])*M_PER_MPC;
    z_min=(input_vals[2]-input_vals[3])*M_PER_MPC;
    z_max=(input_vals[2]+input_vals[3])*M_PER_MPC;
  }

  SID_log("Marking particles in GADGET file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);

  // Read GADGET file
  SID_log("Reading GADGET file...",SID_LOG_OPEN|SID_LOG_TIMER);
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  ADaPS_store(&(plist.data),(void *)(&x_min),"rank_x_min",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),(void *)(&x_max),"rank_x_max",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),(void *)(&y_min),"rank_y_min",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),(void *)(&y_max),"rank_y_max",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),(void *)(&z_min),"rank_z_min",ADaPS_SCALAR_DOUBLE);
  ADaPS_store(&(plist.data),(void *)(&z_max),"rank_z_max",ADaPS_SCALAR_DOUBLE);
  read_gadget_binary(filename_in,snapshot,&plist,READ_GADGET_NO_HUBBLE);
  SID_log("Done.",SID_LOG_CLOSE);

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
