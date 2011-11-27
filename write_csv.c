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
  GBPREAL    dummy_f=0.;
  GBPREAL   *M_gas;
  GBPREAL   *x_gas;
  GBPREAL   *y_gas;
  GBPREAL   *z_gas;
  GBPREAL   *vx_gas;
  GBPREAL   *vy_gas;
  GBPREAL   *vz_gas;
  GBPREAL   *rho_gas;
  GBPREAL   *T_gas;
  GBPREAL   *s_gas;
  GBPREAL   *Z_gas;
  GBPREAL   *phi_gas;
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
  GBPREAL   *M_dark;
  GBPREAL   *x_dark;
  GBPREAL   *y_dark;
  GBPREAL   *z_dark;
  GBPREAL   *vx_dark;
  GBPREAL   *vy_dark;
  GBPREAL   *vz_dark;
  GBPREAL   *s_dark;
  GBPREAL   *rho_dark;
  GBPREAL   *phi_dark;
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
  GBPREAL   *M_star;
  GBPREAL   *x_star;
  GBPREAL   *y_star;
  GBPREAL   *z_star;
  GBPREAL   *vx_star;
  GBPREAL   *vy_star;
  GBPREAL   *vz_star;
  GBPREAL   *Z_star;
  GBPREAL   *t_form_star;
  GBPREAL   *s_star;
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
    M_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"M_gas");
    M_gas_flag=TRUE;
  }
  else
    M_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"x_gas")){
    x_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"x_gas");
    x_gas_flag=TRUE;
  }
  else
    x_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"y_gas")){
    y_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"y_gas");
    y_gas_flag=TRUE;
  }
  else
    y_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"z_gas")){
    z_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"z_gas");
    z_gas_flag=TRUE;
  }
  else
    z_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"vx_gas")){
    vx_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"vx_gas");
    vx_gas_flag=TRUE;
  }
  else
    vx_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"vy_gas")){
    vy_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"vy_gas");
    vy_gas_flag=TRUE;
  }
  else
    vy_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"vz_gas")){
    vz_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"vz_gas");
    vz_gas_flag=TRUE;
  }
  else
    vz_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"rho_gas")){
    rho_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"rho_gas");
    rho_gas_flag=TRUE;
  }
  else
    rho_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"T_gas")){
    T_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"T_gas");
    T_gas_flag=TRUE;
  }
  else
    T_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"s_gas")){
    s_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"s_gas");
    s_gas_flag=TRUE;
  }
  else
    s_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"Z_gas")){
    Z_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"Z_gas");
    Z_gas_flag=TRUE;
  }
  else
    Z_gas_flag=FALSE;
  if(ADaPS_exist(plist->data,"phi_gas")){
    phi_gas     =(GBPREAL *)ADaPS_fetch(plist->data,"phi_gas");
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
    M_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"M_dark");
    M_dark_flag=TRUE;
  }
  else
    M_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"x_dark")){
    x_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"x_dark");
    x_dark_flag=TRUE;
  }
  else
    x_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"y_dark")){
    y_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"y_dark");
    y_dark_flag=TRUE;
  }
  else
    y_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"z_dark")){
    z_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"z_dark");
    z_dark_flag=TRUE;
  }
  else
    z_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"vx_dark")){
    vx_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"vx_dark");
    vx_dark_flag=TRUE;
  }
  else
    vx_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"vy_dark")){
    vy_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"vy_dark");
    vy_dark_flag=TRUE;
  }
  else
    vy_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"vz_dark")){
    vz_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"vz_dark");
    vz_dark_flag=TRUE;
  }
  else
    vz_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"rho_dark")){
    rho_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"s_dark");
    rho_dark_flag=TRUE;
  }
  else
    rho_dark_flag=FALSE;
  if(ADaPS_exist(plist->data,"phi_dark")){
    phi_dark     =(GBPREAL *)ADaPS_fetch(plist->data,"phi_dark");
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
    M_star     =(GBPREAL *)ADaPS_fetch(plist->data,"M_star");
    M_star_flag=TRUE;
  }
  else
    M_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"x_star")){
    x_star     =(GBPREAL *)ADaPS_fetch(plist->data,"x_star");
    x_star_flag=TRUE;
  }
  else
    x_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"y_star")){
    y_star     =(GBPREAL *)ADaPS_fetch(plist->data,"y_star");
    y_star_flag=TRUE;
  }
  else
    y_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"z_star")){
    z_star     =(GBPREAL *)ADaPS_fetch(plist->data,"z_star");
    z_star_flag=TRUE;
  }
  else
    z_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"vx_star")){
    vx_star     =(GBPREAL *)ADaPS_fetch(plist->data,"vx_star");
    vx_star_flag=TRUE;
  }
  else
    vx_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"vy_star")){
    vy_star     =(GBPREAL *)ADaPS_fetch(plist->data,"vy_star");
    vy_star_flag=TRUE;
  }
  else
    vy_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"vz_star")){
    vz_star     =(GBPREAL *)ADaPS_fetch(plist->data,"vz_star");
    vz_star_flag=TRUE;
  }
  else
    vz_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"Z_star")){
    Z_star     =(GBPREAL *)ADaPS_fetch(plist->data,"Z_star");
    Z_star_flag=TRUE;
  }
  else
    Z_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"t_form_star")){
    t_form_star     =(GBPREAL *)ADaPS_fetch(plist->data,"t_form_star");
    t_form_star_flag=TRUE;
  }
  else
    t_form_star_flag=FALSE;
  if(ADaPS_exist(plist->data,"s_star")){
    s_star     =(GBPREAL *)ADaPS_fetch(plist->data,"s_star");
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
  SID_Barrier(SID.COMM_WORLD);
  }

}
