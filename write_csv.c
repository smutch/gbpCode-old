#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void write_csv(char       *filename_out,
               plist_info *plist){
  double  time;
  int     n_dim;
  int     n_gas;
  int     n_dark;
  int     n_star;
  int     n_particles;
  int     i;
  double  h_Hubble;
  REAL    dummy_f=0.;
  REAL   *M_gas;
  REAL   *x_gas;
  REAL   *y_gas;
  REAL   *z_gas;
  REAL   *vx_gas;
  REAL   *vy_gas;
  REAL   *vz_gas;
  REAL   *rho_gas;
  REAL   *T_gas;
  REAL   *s_gas;
  REAL   *Z_gas;
  REAL   *phi_gas;
  int     M_gas_flag;
  int     x_gas_flag;
  int     y_gas_flag;
  int     z_gas_flag;
  int     vx_gas_flag;
  int     vy_gas_flag;
  int     vz_gas_flag;
  int     rho_gas_flag;
  int     T_gas_flag;
  int     s_gas_flag;
  int     Z_gas_flag;
  int     phi_gas_flag;
  REAL   *M_dark;
  REAL   *x_dark;
  REAL   *y_dark;
  REAL   *z_dark;
  REAL   *vx_dark;
  REAL   *vy_dark;
  REAL   *vz_dark;
  REAL   *s_dark;
  REAL   *rho_dark;
  REAL   *phi_dark;
  int     M_dark_flag;
  int     x_dark_flag;
  int     y_dark_flag;
  int     z_dark_flag;
  int     vx_dark_flag;
  int     vy_dark_flag;
  int     vz_dark_flag;
  int     s_dark_flag;
  int     rho_dark_flag;
  int     phi_dark_flag;
  REAL   *M_star;
  REAL   *x_star;
  REAL   *y_star;
  REAL   *z_star;
  REAL   *vx_star;
  REAL   *vy_star;
  REAL   *vz_star;
  REAL   *Z_star;
  REAL   *t_form_star;
  REAL   *s_star;
  int     M_star_flag;
  int     x_star_flag;
  int     y_star_flag;
  int     z_star_flag;
  int     vx_star_flag;
  int     vy_star_flag;
  int     vz_star_flag;
  int     Z_star_flag;
  int     t_form_star_flag;
  int     s_star_flag;
  FILE   *fp;
  int     i_rank;
  int     n_species;

  for(i_rank=0;i_rank<SID.n_proc;i_rank++){
  if(i_rank==SID.My_rank){

  h_Hubble=((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0];

  /*************/
  /* Open file */
  /*************/
  if(i_rank==MASTER_RANK)
    fp=fopen(filename_out,"w");
  else
    fp=fopen(filename_out,"a");

  /****************/
  /* Write header */
  /****************/
  /* Time */
  if(ADaPS_exist(plist->data,"time"))
    time=((double *)ADaPS_fetch(plist->data,"time"))[0];
  else if(ADaPS_exist(plist->data,"expansion_factor"))
    time=((double *)ADaPS_fetch(plist->data,"expansion_factor"))[0];
  else if(ADaPS_exist(plist->data,"redshift"))
    time=((double *)ADaPS_fetch(plist->data,"redshift"))[0];
  else
    time=0.0;

  /* No. of dimensions and particles */
  if(ADaPS_exist(plist->data,"n_dims"))
    n_dim=((int *)ADaPS_fetch(plist->data,"n_dims"))[0];
  else
    n_dim=3;
  n_species=0;
  if(ADaPS_exist(plist->data,"n_gas"))
    n_gas=((int *)ADaPS_fetch(plist->data,"n_gas"))[0];
  else
    n_gas=0;
  if(n_gas>0)  n_species++;
  if(ADaPS_exist(plist->data,"n_dark"))
    n_dark=((int *)ADaPS_fetch(plist->data,"n_dark"))[0];
  else
    n_dark=0;
  if(n_dark>0) n_species++;
  if(ADaPS_exist(plist->data,"n_star"))
    n_star=((int *)ADaPS_fetch(plist->data,"n_star"))[0];
  else
    n_star=0;
  if(n_star>0) n_species++;
  n_particles=n_gas+n_dark+n_star;

  // csv header
  fprintf(fp,"float32 Position[0], float32 Position[1], float32 Position[2], float16 Velocity[0], float16 Velocity[1], float16 Velocity[2]\n");

  /***********************/
  /* Write gas particles */
  /***********************/
  if(ADaPS_exist(plist->data,"M_gas")){
    M_gas     =(REAL *)ADaPS_fetch(plist->data,"M_gas");
    M_gas_flag=TRUE;
  }
  else
    M_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"x_gas")){
    x_gas     =(REAL *)ADaPS_fetch(plist->data,"x_gas");
    x_gas_flag=TRUE;
  }
  else
    x_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"y_gas")){
    y_gas     =(REAL *)ADaPS_fetch(plist->data,"y_gas");
    y_gas_flag=TRUE;
  }
  else
    y_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"z_gas")){
    z_gas     =(REAL *)ADaPS_fetch(plist->data,"z_gas");
    z_gas_flag=TRUE;
  }
  else
    z_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"vx_gas")){
    vx_gas     =(REAL *)ADaPS_fetch(plist->data,"vx_gas");
    vx_gas_flag=TRUE;
  }
  else
    vx_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"vy_gas")){
    vy_gas     =(REAL *)ADaPS_fetch(plist->data,"vy_gas");
    vy_gas_flag=TRUE;
  }
  else
    vy_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"vz_gas")){
    vz_gas     =(REAL *)ADaPS_fetch(plist->data,"vz_gas");
    vz_gas_flag=TRUE;
  }
  else
    vz_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"rho_gas")){
    rho_gas     =(REAL *)ADaPS_fetch(plist->data,"rho_gas");
    rho_gas_flag=TRUE;
  }
  else
    rho_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"T_gas")){
    T_gas     =(REAL *)ADaPS_fetch(plist->data,"T_gas");
    T_gas_flag=TRUE;
  }
  else
    T_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"s_gas")){
    s_gas     =(REAL *)ADaPS_fetch(plist->data,"s_gas");
    s_gas_flag=TRUE;
  }
  else
    s_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"Z_gas")){
    Z_gas     =(REAL *)ADaPS_fetch(plist->data,"Z_gas");
    Z_gas_flag=TRUE;
  }
  else
    Z_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"phi_gas")){
    phi_gas     =(REAL *)ADaPS_fetch(plist->data,"phi_gas");
    phi_gas_flag=TRUE;
  }
  else
    phi_gas_flag=FALSE;
  for(i=0;i<n_gas;i++){
    if(n_species>1) fprintf(fp,"0");
    // x 
    if(x_gas_flag) dummy_f=(float)(x_gas[i]*h_Hubble/M_PER_MPC);
    else           dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // y 
    if(y_gas_flag) dummy_f=(float)(y_gas[i]*h_Hubble/M_PER_MPC);
    else           dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // z 
    if(z_gas_flag) dummy_f=(float)(z_gas[i]*h_Hubble/M_PER_MPC);
    else           dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // vx 
    if(vx_gas_flag) dummy_f=(float)(vx_gas[i]/1e3);
    else            dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // vy 
    if(vy_gas_flag) dummy_f=(float)(vy_gas[i]/1e3);
    else            dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // vz 
    if(vz_gas_flag) dummy_f=(float)(vz_gas[i]/1e3);
    else            dummy_f=0.;
/*
    fprintf(fp," %11.4le",dummy_f);
    // rho 
    if(rho_gas_flag) dummy_f=(float)(rho_gas[i]);
    else             dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // T 
    if(rho_gas_flag) dummy_f=(float)(T_gas[i]);
    else             dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
*/
    fprintf(fp,"\n");
  }

  /***********************/
  /* Write dark particles */
  /***********************/
  if(ADaPS_exist(plist->data,"M_dark")){
    M_dark     =(REAL *)ADaPS_fetch(plist->data,"M_dark");
    M_dark_flag=TRUE;
  }
  else
    M_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"x_dark")){
    x_dark     =(REAL *)ADaPS_fetch(plist->data,"x_dark");
    x_dark_flag=TRUE;
  }
  else
    x_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"y_dark")){
    y_dark     =(REAL *)ADaPS_fetch(plist->data,"y_dark");
    y_dark_flag=TRUE;
  }
  else
    y_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"z_dark")){
    z_dark     =(REAL *)ADaPS_fetch(plist->data,"z_dark");
    z_dark_flag=TRUE;
  }
  else
    z_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"vx_dark")){
    vx_dark     =(REAL *)ADaPS_fetch(plist->data,"vx_dark");
    vx_dark_flag=TRUE;
  }
  else
    vx_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"vy_dark")){
    vy_dark     =(REAL *)ADaPS_fetch(plist->data,"vy_dark");
    vy_dark_flag=TRUE;
  }
  else
    vy_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"vz_dark")){
    vz_dark     =(REAL *)ADaPS_fetch(plist->data,"vz_dark");
    vz_dark_flag=TRUE;
  }
  else
    vz_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"rho_dark")){
    rho_dark     =(REAL *)ADaPS_fetch(plist->data,"s_dark");
    rho_dark_flag=TRUE;
  }
  else
    rho_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"phi_dark")){
    phi_dark     =(REAL *)ADaPS_fetch(plist->data,"phi_dark");
    phi_dark_flag=TRUE;
  }
  else
    phi_dark_flag=FALSE;
  for(i=0;i<n_dark;i++){
    if(n_species>1) fprintf(fp,"1");
    // x 
    if(x_dark_flag) dummy_f=(float)(x_dark[i]*h_Hubble/M_PER_MPC);
    else            dummy_f=0.;
    //fprintf(fp," %11.4le,",dummy_f);
    fprintf(fp," %5.3lf,",dummy_f);
    // y 
    if(y_dark_flag) dummy_f=(float)(y_dark[i]*h_Hubble/M_PER_MPC);
    else            dummy_f=0.;
    //fprintf(fp," %11.4le,",dummy_f);
    fprintf(fp," %5.3lf,",dummy_f);
    // z 
    if(z_dark_flag) dummy_f=(float)(z_dark[i]*h_Hubble/M_PER_MPC);
    else            dummy_f=0.;
    //fprintf(fp," %11.4le,",dummy_f);
    fprintf(fp," %5.3lf,",dummy_f);
    // vx 
    if(vx_dark_flag) dummy_f=(float)(vx_dark[i]/1e3);
    else             dummy_f=0.;
    //fprintf(fp," %11.4le,",dummy_f);
    fprintf(fp," %5.3lf,",dummy_f);
    // vy
    if(vy_dark_flag) dummy_f=(float)(vy_dark[i]/1e3);
    else             dummy_f=0.;
    //fprintf(fp," %11.4le,",dummy_f);
    fprintf(fp," %5.3lf,",dummy_f);
    // vz 
    if(vz_dark_flag) dummy_f=(float)(vz_dark[i]/1e3);
    else             dummy_f=0.;
    //fprintf(fp," %11.4le",dummy_f);
    fprintf(fp," %5.3lf",dummy_f);
/*
    // rho 
    if(rho_dark_flag) dummy_f=(float)(rho_dark[i]);
    else              dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // phi 
    if(phi_dark_flag) dummy_f=(float)(phi_dark[i]);
    else              dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
*/
    fprintf(fp,"\n");
  }

  /***********************/
  /* Write star particles */
  /***********************/
  if(ADaPS_exist(plist->data,"M_star")){
    M_star     =(REAL *)ADaPS_fetch(plist->data,"M_star");
    M_star_flag=TRUE;
  }
  else
    M_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"x_star")){
    x_star     =(REAL *)ADaPS_fetch(plist->data,"x_star");
    x_star_flag=TRUE;
  }
  else
    x_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"y_star")){
    y_star     =(REAL *)ADaPS_fetch(plist->data,"y_star");
    y_star_flag=TRUE;
  }
  else
    y_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"z_star")){
    z_star     =(REAL *)ADaPS_fetch(plist->data,"z_star");
    z_star_flag=TRUE;
  }
  else
    z_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"vx_star")){
    vx_star     =(REAL *)ADaPS_fetch(plist->data,"vx_star");
    vx_star_flag=TRUE;
  }
  else
    vx_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"vy_star")){
    vy_star     =(REAL *)ADaPS_fetch(plist->data,"vy_star");
    vy_star_flag=TRUE;
  }
  else
    vy_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"vz_star")){
    vz_star     =(REAL *)ADaPS_fetch(plist->data,"vz_star");
    vz_star_flag=TRUE;
  }
  else
    vz_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"Z_star")){
    Z_star     =(REAL *)ADaPS_fetch(plist->data,"Z_star");
    Z_star_flag=TRUE;
  }
  else
    Z_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"t_form_star")){
    t_form_star     =(REAL *)ADaPS_fetch(plist->data,"t_form_star");
    t_form_star_flag=TRUE;
  }
  else
    t_form_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"s_star")){
    s_star     =(REAL *)ADaPS_fetch(plist->data,"s_star");
    s_star_flag=TRUE;
  }
  else
    s_star_flag=FALSE;
  for(i=0;i<n_star;i++){
    if(n_species>1) fprintf(fp,"2");
    // x 
    if(x_star_flag) dummy_f=(float)(x_star[i]*h_Hubble/M_PER_MPC);
    else            dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // y 
    if(y_star_flag) dummy_f=(float)(y_star[i]*h_Hubble/M_PER_MPC);
    else            dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // z 
    if(z_star_flag) dummy_f=(float)(z_star[i]*h_Hubble/M_PER_MPC);
    else            dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // vx 
    if(vx_star_flag) dummy_f=(float)(vx_star[i]/1e3);
    else             dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // vy 
    if(vy_star_flag) dummy_f=(float)(vy_star[i]/1e3);
    else             dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // vz 
    if(vz_star_flag) dummy_f=(float)(vz_star[i]/1e3);
    else             dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
/*
    // Z 
    if(Z_star_flag) dummy_f=(float)(Z_star[i]);
    else            dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
    // t_form 
    if(t_form_star_flag) dummy_f=(float)(t_form_star[i]);
    else                 dummy_f=0.;
    fprintf(fp," %11.4le",dummy_f);
*/
    fprintf(fp,"\n");
  }

  /******************/
  /* Close the file */
  /******************/
  fclose(fp);
  }
  #ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  }

}
