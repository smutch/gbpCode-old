#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo_bias.h>
double bias_model_BPR_integral(cosmo_info **cosmo,
                               double       z);
double bias_model_BPR_integral(cosmo_info **cosmo,
                               double       z){
   double       z_max=10000.;
   interp_info *interp;
   if(!ADaPS_exist(*cosmo,"bias_model_BPR_Iz_interp")){
      int     n_int;
      int     i_int;
      double  dz;
      double  Omega_M,Omega_k,Omega_Lambda;
      double  z_temp;
      double *x_int;
      double *y_int;
      double  log_z;
      double  dlog_z;
      n_int       =250;
      Omega_M     =((double *)ADaPS_fetch(*cosmo,"Omega_M"))[0];
      Omega_k     =((double *)ADaPS_fetch(*cosmo,"Omega_k"))[0];
      Omega_Lambda=((double *)ADaPS_fetch(*cosmo,"Omega_Lambda"))[0];
      x_int       = (double *)SID_malloc(sizeof(double)*n_int);
      y_int       = (double *)SID_malloc(sizeof(double)*n_int);
      i_int=0;
      x_int[i_int]=0.;
      y_int[i_int]=pow((1.+x_int[i_int])/E_z(Omega_M,Omega_k,Omega_Lambda,x_int[i_int]),3.);
      i_int++;
      x_int[i_int]=take_log10(z_max)/(double)(n_int-1);
      y_int[i_int]=pow((1.+x_int[i_int])/E_z(Omega_M,Omega_k,Omega_Lambda,x_int[i_int]),3.);
      log_z    =take_log10(x_int[i_int]);
      dlog_z   =(take_log10(z_max)-log_z)/(double)(n_int-2);
      for(i_int++,log_z+=dlog_z;i_int<(n_int-1);i_int++,log_z+=dlog_z){
        x_int[i_int]=take_alog10(log_z);
        y_int[i_int]=pow((1.+x_int[i_int])/E_z(Omega_M,Omega_k,Omega_Lambda,x_int[i_int]),3.);
      }
      x_int[i_int]=z_max;
      y_int[i_int]=pow((1.+x_int[i_int])/E_z(Omega_M,Omega_k,Omega_Lambda,x_int[i_int]),3.);
      init_interpolate(x_int,y_int,(size_t)n_int,gsl_interp_cspline,&interp);
      SID_free(SID_FARG x_int);
      SID_free(SID_FARG y_int);
      ADaPS_store_interp(cosmo,
                         (void *)(interp),
                         "bias_model_BPR_Iz_interp");
          
   }
   else
      interp=(interp_info *)ADaPS_fetch(*cosmo,"bias_model_BPR_Iz_interp");
   return(interpolate_integral(interp,z,z_max));
}

double bias_model(double       x_in,
                  double       delta_c,
                  double       z,
                  cosmo_info **cosmo,
                  int          mode){

   // Decide what the input is
   int    flag_Vmax_ordinate=check_mode_for_flag(mode,BIAS_MODEL_VMAX_ORDINATE);
   double M_R;
   double V_max;
   if(flag_Vmax_ordinate)
      V_max=x_in;
   else
      M_R=x_in;

   double bias;
   int    flag_done=FALSE;

   // Tinker et al 2010
   if(check_mode_for_flag(mode,BIAS_MODEL_TRK)){
      if(flag_done)
         SID_trap_error("Mode flag (%d) is invalid in bias_model().  Multiple model definitions.",ERROR_LOGIC,mode);
      flag_done=TRUE;
      if(flag_Vmax_ordinate)
         M_R=Vmax_to_Mvir_NFW(V_max,z,NFW_MODE_DEFAULT,cosmo);
      double y    =take_log10(Delta_vir(z,*cosmo));
      double A    =1.+0.24*y*exp(-pow(4./y,4.));
      double B    =0.183;
      double C    =0.019+0.107*y+0.19*exp(-pow(4./y,4.));
      double a    =0.44*y-0.88;
      double b    =1.5;
      double c    =2.4;
      double sigma=sigma_M(cosmo,M_R,z,PSPEC_LINEAR_TF,PSPEC_ALL_MATTER);
      double nu   =delta_c/sigma;
      bias=1.-A*pow(nu,a)/(pow(nu,a)+pow(delta_c,a))+B*pow(nu,b)+C*pow(nu,c);
   }
   // Basilakos and Plionis (2001;2003) w/ Papageorgiou et al 2013 coeeficients
   if(check_mode_for_flag(mode,BIAS_MODEL_BPR)){
      if(flag_done)
         SID_trap_error("Mode flag (%d) is invalid in bias_model().  Multiple model definitions.",ERROR_LOGIC,mode);
      flag_done=TRUE;
      if(flag_Vmax_ordinate)
         M_R=Vmax_to_Mvir_NFW(V_max,z,NFW_MODE_DEFAULT,cosmo);
      double alpha_1= 4.53;
      double alpha_2=-0.41;
      double beta_1 = 0.37;
      double beta_2 = 0.36;
      double I_z;
      double C_1;
      double C_2;
      double Omega_M,Omega_k,Omega_Lambda,h_Hubble;
      double Ez;
      Omega_M     =((double *)ADaPS_fetch(*cosmo,"Omega_M"))[0];
      Omega_k     =((double *)ADaPS_fetch(*cosmo,"Omega_k"))[0];
      Omega_Lambda=((double *)ADaPS_fetch(*cosmo,"Omega_Lambda"))[0];
      h_Hubble    =((double *)ADaPS_fetch(*cosmo,"h_Hubble"))[0];
      Ez          =E_z(Omega_M,Omega_k,Omega_Lambda,z);
      I_z         =bias_model_BPR_integral(cosmo,z);
      C_1         =alpha_1*pow(M_R/(1e13*M_SOL/h_Hubble),beta_1);
      C_2         =alpha_2*pow(M_R/(1e13*M_SOL/h_Hubble),beta_2);
      bias=(C_1+C_2*I_z)*Ez+1.;
   }
   // Poole et al 2014
   if(check_mode_for_flag(mode,BIAS_MODEL_POOLE_HALO)){
      if(flag_done)
         SID_trap_error("Mode flag (%d) is invalid in bias_model().  Multiple model definitions.",ERROR_LOGIC,mode);
      flag_done=TRUE;
      if(!flag_Vmax_ordinate)
         V_max = V_max_NFW(M_R,z,NFW_MODE_DEFAULT,cosmo)*1e-3;
      else
         V_max*=1e-3;
      double V_SF_o;
      double V_SF_z;
      double s_V_o ;
      double s_V_z ;
      double b_o_o ;
      double b_o_z ;
      double b_V_o ;
      double b_V_z ;
      if(check_mode_for_flag(mode,BIAS_MODEL_POOLE_SUBSTRUCTURE)){
         V_SF_o= 5.326176e-02;
         V_SF_z=-1.673868e-01;
         s_V_o = 4.026941e-01;
         s_V_z = 6.096567e-01;
         b_o_o =-1.974311e-01;
         b_o_z = 2.138219e-01;
         b_V_o = 2.707540e-01;
         b_V_z = 8.202001e-02;
      }
      else{
         V_SF_o= 2.819063e-02;
         V_SF_z=-1.381993e-01;
         s_V_o = 3.685953e-01;
         s_V_z = 6.154695e-01;
         b_o_o =-3.793559e-01;
         b_o_z = 3.074326e-01;
         b_V_o = 3.147507e-01;
         b_V_z = 6.072666e-02;
      }
      double b_o   = b_o_o+b_o_z*z;
      double b_V   = (b_V_o+b_V_z*z)/220.;
      double V_SF  = take_alog10(V_SF_o+V_SF_z*z)*220.;
      double s_V   = (s_V_o+s_V_z*z)/220.;
      double s     = s_V*fabs(V_max-V_SF);
      bias         = take_alog10(0.5*(b_o+b_V*V_max));
   }
   if(check_mode_for_flag(mode,BIAS_MODEL_POOLE_ZSPACE)){
      if(flag_done)
         SID_trap_error("Mode flag (%d) is invalid in bias_model().  Multiple model definitions.",ERROR_LOGIC,mode);
      flag_done=TRUE;
      if(!flag_Vmax_ordinate)
         V_max = V_max_NFW(M_R,z,NFW_MODE_DEFAULT,cosmo)*1e-3;
      else
         V_max*=1e-3;
      double V_SF_o;
      double V_SF_z;
      double s_V_o ;
      double s_V_zz;
      double b_o_o ;
      double b_o_zz;
      double b_V_o ;
      double b_V_z ;
      double z_b_c ;
      if(check_mode_for_flag(mode,BIAS_MODEL_POOLE_SUBSTRUCTURE)){
         V_SF_o= 3.173152e-01;
         V_SF_z=-1.599133e-01;
         s_V_o = 5.344408e-01;
         s_V_zz= 7.102406e-02;
         b_o_o = 2.198795e-01;
         b_o_zz=-3.749491e-02;
         b_V_o =-4.628602e-02;
         b_V_z =-1.832620e-02;
         z_b_c = 9.292014e-01;
      }
      else{
         V_SF_o= 3.100167e-01;
         V_SF_z=-2.026411e-01;
         s_V_o = 3.342258e-01;
         s_V_zz= 9.233431e-02;
         b_o_o = 2.206163e-01;
         b_o_zz=-4.419126e-02;
         b_V_o =-4.804747e-02;
         b_V_z =-1.454479e-02;
         z_b_c = 7.852721e-01;
      }
      double b_o   = b_o_o+b_o_zz*(z-z_b_c)*(z-z_b_c);
      double b_V   = (b_V_o+b_V_z*z)/220.;
      double V_SF  = 220.*take_alog10(V_SF_o+V_SF_z*z);
      double s_V   = (s_V_o+s_V_zz*z*z)/220.;
      double s     = s_V*fabs(V_max-V_SF);
      bias         = take_alog10(0.5*(b_o+b_V*V_max));
   }
   if(check_mode_for_flag(mode,BIAS_MODEL_POOLE_TOTAL)){
      if(flag_done)
         SID_trap_error("Mode flag (%d) is invalid in bias_model().  Multiple model definitions.",ERROR_LOGIC,mode);
      flag_done=TRUE;
      if(!flag_Vmax_ordinate)
         V_max = V_max_NFW(M_R,z,NFW_MODE_DEFAULT,cosmo)*1e-3;
      else
         V_max*=1e-3;
      double V_SF_o;
      double V_SF_z;
      double s_V_o ;
      double s_V_z ;
      double b_o_o ;
      double b_o_z ;
      double b_V_o ;
      double b_V_z ;
      if(check_mode_for_flag(mode,BIAS_MODEL_POOLE_SUBSTRUCTURE)){
         V_SF_o= 2.128681e-01;
         V_SF_z=-2.280569e-01;
         s_V_o = 1.118792e+00;
         s_V_z = 6.121383e-01;
         b_o_o = 1.198136e-02;
         b_o_z = 2.112670e-01;
         b_V_o = 2.153513e-01;
         b_V_z = 7.763461e-02;
      }
      else{
         V_SF_o= 2.041659e-01;
         V_SF_z=-2.966696e-01;
         s_V_o = 9.408169e-01;
         s_V_z = 4.514711e-01;
         b_o_o =-1.534952e-01;
         b_o_z = 2.799483e-01;
         b_V_o = 2.547096e-01;
         b_V_z = 6.760491e-02;
      }
      double b_o   = b_o_o+b_o_z*z;
      double b_V   = (b_V_o+b_V_z*z)/220.;
      double V_SF  = V_SF_o+V_SF_z*z;
      double s_V   = (s_V_o+s_V_z*z)/220.;
      double s     = s_V*fabs(V_max-V_SF);
      bias         = take_alog10(0.5*(b_o+b_V*V_max));
   }
   if(!flag_done)
      SID_trap_error("Mode flag (%d) is invalid in bias_model().  No model definition.",ERROR_LOGIC,mode);
   // Apply the Kaiser '87 model to whatever model has been processed above.  Be careful, there
   //   are some mode flags (such as BIAS_MODE_POOLE_ZSPACE) for which this does not make sence 
   //   and we presently don't check for this.
   if(check_mode_for_flag(mode,BIAS_MODEL_KAISER_BOOST)){
      double Omega_M_z=Omega_z(z,(*cosmo));
      double f        =pow(Omega_M_z,0.55);
      double beta     =f/bias;
      double boost    =pow(1.+TWO_THIRDS*beta+0.2*beta*beta,0.5);
      bias=boost;
   }
   if(check_mode_for_flag(mode,BIAS_MODEL_KAISER)){
      double Omega_M_z=Omega_z(z,(*cosmo));
      double f        =pow(Omega_M_z,0.55);
      double beta     =f/bias;
      double boost    =pow(1.+TWO_THIRDS*beta+0.2*beta*beta,0.5);
      bias*=boost;
   }
   return(bias);
}

