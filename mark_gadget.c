#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  int         i;
  plist_info  plist;
  char        filename_in[256];
  char        filename_out[256];
  int         flag_mode;
  int         n_input_vals;
  size_t      n_marked;
  double     *input_vals;
  double      h_Hubble;

  SID_init(&argc,&argv,NULL);

  // Parse command line
  if(argc<=1){
    SID_log_error("\nsyntax: %s filename_in filename_out [box/sphere] [x_min,xmax,y_min,y_max,z_min,z_max/x_cen,y_cen,z_cen,radius]",ERROR_SYNTAX,argv[0]);
    SID_trap_error("------\n",ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in, argv[1]);
    strcpy(filename_out,argv[2]);
    if(!strcmp(argv[3],"box"))
      flag_mode=VOLUME_BOX;
    else if(!strcmp(argv[3],"sphere"))
      flag_mode=VOLUME_SPHERE;
    else
      SID_trap_error("Valid modes are \"box\" or \"sphere\".",ERROR_SYNTAX);
    if(flag_mode==VOLUME_BOX)
      n_input_vals=6;
    else if(flag_mode=VOLUME_SPHERE)
      n_input_vals=4;
  }

  SID_log("Marking particles in GADGET file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in);

  // Read GADGET file
  SID_log("Reading GADGET file...",SID_LOG_OPEN|SID_LOG_TIMER);
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_gadget_binary(filename_in,&plist,READ_GADGET_DEFAULT);
  h_Hubble=((double *)ADaPS_fetch(plist.data,"h_Hubble"))[0];
  SID_log("Done.",SID_LOG_CLOSE);

  // Initialize input parameters
  input_vals=(double *)SID_malloc(sizeof(double)*n_input_vals);
  for(i=0;i<n_input_vals;i++)
    input_vals[i]=(double)atof(argv[4+i])*M_PER_MPC/h_Hubble;

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
