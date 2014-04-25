/**********************/
/* Cosmology routines */
/**********************/
#ifndef COSMO_AWAKE
#define COSMO_AWAKE

#define PSPEC_LINEAR_TF       0
#define PSPEC_LINEAR_BBKS     1
#define PSPEC_NONLINEAR_JAIN  2
#define PSPEC_NONLINEAR_PD    3
#define PSPEC_NONLINEAR_SMITH 4

#define PSPEC_ALL_MATTER  0
#define PSPEC_DARK_MATTER 1
#define PSPEC_BARYON      2

#define MF_JENKINS        0
#define MF_PS             1
#define MF_ST             2

#define DELTA_BRYAN_NORMAN 0
#define DELTA_LINEAR       1

#define DELTAT_A_MIN_A     1e-5

typedef ADaPS cosmo_info;

double a_of_z(double z);
double z_of_a(double a);
void init_deltat_a(cosmo_info **cosmo);
double deltat_a(cosmo_info **cosmo,double a_1,double a_2);
double R_Delta_z(double      M_Delta,
		 double      Delta,
		 double      redshift,
                 cosmo_info *cosmo);
double Omega_z(double     redshift,cosmo_info *cosmo);
double Delta_vir(double redshift,cosmo_info *cosmo);
double Ha_Ho(double a,cosmo_info *cosmo);
double H_z(double redshift,cosmo_info *cosmo);
double H_convert(double Hz);
double dlna_dtau(double      a,
		 cosmo_info *cosmo);
double rho_crit_z(double redshift, cosmo_info *cosmo);
double rho_crit_z_strip(double redshift,
                        double h_Hubble,
                        double Omega_M, 
                        double Omega_Lambda);
double D_Hubble(double h_Hubble);
double E_z(double Omega_M, double Omega_k, double Omega_Lambda, double z);
double D_comove(double z,cosmo_info *cosmo);
double D_comove_transverse(double z,cosmo_info *cosmo);
double D_angular(double z,cosmo_info *cosmo);
double D_angular_1to2(double z_1,double z_2,cosmo_info *cosmo);
double D_luminosity(double z,cosmo_info *cosmo);
void init_cosmo_std(cosmo_info **cosmo);
void init_cosmo(cosmo_info **cosmo,
		  double       Omega_Lambda,
		  double       Omega_M,
		  double       Omega_k,
                  double       Omega_b,
                  double       f_gas,
		  double       h_Hubble,
		  double       sigma_8,
		  double       n_spectral);
void free_cosmo(cosmo_info **cosmo);

// Bias model stuff
#define BIAS_MODEL_TRK                  2
#define BIAS_MODEL_BPR                  4
#define BIAS_MODEL_POOLE_SUBSTRUCTURE   8
#define BIAS_MODEL_POOLE_HALO          16
#define BIAS_MODEL_POOLE_ZSPACE        32
#define BIAS_MODEL_POOLE_TOTAL         64
#define BIAS_MODEL_VMAX_ORDINATE      128
#define BIAS_MODEL_KAISER_BOOST       256
#define BIAS_MODEL_KAISER             512

double bias_model(double       x_in,
                  double       delta_c,
                  double       z,
                  cosmo_info **cosmo,
                  int          mode);

/**************************/
/* Linear theory routines */
/**************************/
double R_of_k(double k);
double k_of_R(double R);
double M_of_R(double R,double z,cosmo_info *cosmo);
double M_of_k(double k,double z,cosmo_info *cosmo);
double R_of_M(double M,double z,cosmo_info *cosmo);
double k_of_M(double M,double z,cosmo_info *cosmo);

void pspec_names(int   mode,
		 int   component,
		 char *mode_name,
		 char *component_name);
void   init_transfer_function(cosmo_info **cosmo);
void   init_power_spectrum_TF(cosmo_info **cosmo);
double power_spectrum(double k, double z, cosmo_info **cosmo, int mode,int component);
void   init_power_spectrum_variance(cosmo_info **cosmo,double z,int mode,int component);
double power_spectrum_normalization(cosmo_info *cosmo,
				    int         mode,
				    int         component);
double power_spectrum_variance(double       k_interp,
			       double       redshift,
			       cosmo_info **cosmo,
			       int          mode,
			       int          component);
double dln_Inv_sigma_dlogM(cosmo_info *cosmo,
			   double      M_interp,
			   double      z,
			   int         mode,
			   int         component);
void init_sigma_M(cosmo_info **cosmo,
                  double       z,
                  int          mode,
                  int          component);
double sigma_M(cosmo_info *cosmo,
               double      M_interp,
               double      z,
               int         mode,
               int         component);
double ln_sigma_M(cosmo_info *cosmo,
                  double      M_interp,
                  double      z,
                  int         mode,
                  int         component);
double ln_Inv_sigma_M(cosmo_info *cosmo,
                  double      M_interp,
                  double      z,
                  int         mode,
                  int         component);
double dln_sigma_dlnM(cosmo_info *cosmo,
		      double      M_interp,
		      double      z,
		      int         mode,
		      int         component);
void evolve_power_spectrum(double       a_start,
                           cosmo_info **cosmo);
double dlogDplus_dloga(double      a,
                       cosmo_info *cosmo);
double dDplus_da(double      a,
                 cosmo_info *cosmo);
double Dplus(double      a,
             cosmo_info *cosmo);
double linear_growth_factor(double       redshift,
                            cosmo_info  *cosmo);
double sigma2_linear(double      M,
                     double      redshift,
                     cosmo_info *cosmo);
double W_k_tophat(double kR);
double M_sc(double       z,
	    cosmo_info **cosmo,
	    int          mode,
	    int          component);
double lk_sc(double       z,
	     cosmo_info **cosmo,
	     int          mode,
	     int          component);

// Mass functions
double scaled_mass_function(double sigma,int select_flag);
double mass_function(double        M_interp,
                     double        z,
                     cosmo_info  **cosmo,
                     int           mode,
                     int           component,
                     int           select_flag);
double mass_function_cumulative(double       M_interp,
                                double       z,
                                cosmo_info  *cosmo,
                                int          mode,
                                int          component,
                                int          select_flag);

/********************************/
/* Nonlinear theory routines    */
/* (Mostly taken from smith2.h) */
/********************************/
double P_NL(double k_NL, double z, cosmo_info *cosmo, int nonlinear);
double P_kappa(double s,cosmo_info *cosmo,int nonlinear);

double P_L_Smith(double a, 
		 double k, 
		 cosmo_info *cosmo, 
		 int    mode);
double Tsqr(double      k,
	    cosmo_info *cosmo);
double int_for_sigma_8(double      k,
		       cosmo_info *cosmo,
		       int         junk1,
		       double      junk2);
double sigma_8_sqr(cosmo_info *cosmo);

double n_L(double a, double k,cosmo_info *cosmo,int nonlinear);
double f_NL(double x, double a, double k, cosmo_info *cosmo, int nonlinear);
double Delta_L_BE(double k, cosmo_info *cosmo, int nonlinear);
double int_for_wint_knl(double logk, cosmo_info *cosmo, int nonlinear,double rglob);
double int_for_wint_neff(double logk, cosmo_info *cosmo, int nonlinear,double rglob);
double int_for_wint_ncur(double logk, cosmo_info *cosmo, int nonlinear,double rglob);
void   wint(double  r, 
	    double *sig, 
	    double *d1, 
	    double *d2, 
	    double  amp, 
	    int     onlysig, 
	    cosmo_info *cosmo, 
	    int     nonlinear);
double slope_NL(double rn, double rncur, double om_m, double om_v);
void   halofit(double rk, 
	       double rn, 
	       double rncur, 
	       double rknl, 
	       double plin, 
	       double om_m, 
	       double om_v, 
	       double *pnl);
double dlog(double x);
double int_for_w(double a,cosmo_info *cosmo,int nonlinear,double glob);
double w(double a, cosmo_info *cosmo, int nonlinear);
double f_K(double w, cosmo_info *cosmo, int nonlinear);
double prob_smith(double z, cosmo_info *cosmo, int nonlinear);
double int_for_g(double aprime,cosmo_info *cosmo,int nonlinear,double aglob);
double g_source(double a, cosmo_info *cosmo, int nonlinear);
double int_for_p_2(double a,cosmo_info *cosmo,int nonlinear,double sglob);



void Omega_a(double  a, 
	     double *omega_m, 
	     double *omega_v, 
	     cosmo_info *cosmo, 
	     int     nonlinear);
double g_smith(double a,cosmo_info *cosmo,int nonlinear);


void error_smith(const char *s);
double dfridr_smith(double  (*func)(double,
				    double,
				    cosmo_info *,
				    int), 
		    double  x, 
		    double  h, 
		    double *err, 
		    double  aa, 
		    cosmo_info *cosmo, 
		    int nonlinear);
void polint_smith(double  xa[], 
		  double  ya[], 
		  int     n, 
		  double  x, 
		  double *y, 
		  double *dy);
double trapzd_smith(double (*func)(double,
				   cosmo_info *,
				   int,
				   double), 
		    double a, 
		    double b, 
		    int    n, 
		    cosmo_info *cosmo, 
		    int    nonlinear,
		    double glob);
double interpol2d(double **f,
		  int nx, double ax, double bx, double dx, double x,
		  int ny, double ay, double by, double dy, double y,
		  double lower, double upper);
double gammln(double xx);
double qromb_smith(double (*func)(double,
				  cosmo_info *,
				  int,
				  double), 
		   double a, 
		   double b,
		   cosmo_info *cosmo,
		   int    nonlinear,
		   double glob);
double qromo_smith(double (*func)(double,
				  cosmo_info *,
				  int,
				  double), 
		   double a, 
		   double b,
		   double (*choose)(double (*)(double,
					       cosmo_info *,
					       int,double), 
				    double, 
				    double, 
				    int,
				    cosmo_info *,
				    int,
				    double),
		   cosmo_info *cosmo,
		   int nonlinear,
		   double glob);
double midpnt_smith(double (*func)(double,
				   cosmo_info *,
				   int,
				   double), 
		    double a, 
		    double b, 
		    int    n, 
		    cosmo_info *cosmo, 
		    int    nonlinear,
		    double glob);

/************************/
/* Routines for NFW etc */
/************************/

#define NFW_MODE_DEFAULT 0

void set_NFW_params(double       M,
		    double       z,
		    int          mode,
		    cosmo_info **cosmo,
		    double      *c_vir,
		    double      *R_vir);
double R_vir_NFW(double       M_vir,
		 double       z,
		 int          mode,
		 cosmo_info **cosmo);
double c_vir_NFW(double       M_vir,
		 double       z,
		 int          mode,
		 cosmo_info **cosmo);
double rho_NFW(double       r,
	       double       M_vir,
	       double       z,
	       int          mode,
	       cosmo_info **cosmo);
double rho_NFW_fft(double       k,
		   double       M_vir,
		   double       z,
		   int          mode,
		   cosmo_info **cosmo);
double M_r_NFW(double       r,
	       double       M_vir,
	       double       z,
	       int          mode,
	       cosmo_info **cosmo);
double V_circ_NFW(double       r,
		  double       M_vir,
		  double       z,
		  int          mode,
		  cosmo_info **cosmo);
double V_circ_vir_NFW(double       M_vir,
		      double       z,
		      int          mode,
		      cosmo_info **cosmo);
double V_max_NFW(double       M_vir,
		 double       z,
		 int          mode,
		 cosmo_info **cosmo);
void init_Vmax_to_Mvir_NFW(cosmo_info **cosmo,
                           int          mode,
                           double       z);
double Vmax_to_Mvir_NFW(double       V_max,
                        double       z,
                        int          mode,
                        cosmo_info **cosmo);
double R_half_V_max_NFW(double       M_vir,
			double       z,
			int          mode,
			cosmo_info **cosmo);
double R_V_max_NFW(double       M_vir,
                   double       z,
                   int          mode,
                   cosmo_info **cosmo);
double Delta_half_V_max_NFW(double       M_vir,
			    double       z,
			    int          mode,
			    cosmo_info **cosmo);
double Delta_half_V_max_NFW(double       M_vir,
			    double       z,
			    int          mode,
			    cosmo_info **cosmo);
double V2_circ(double M, 
	       double r);

/*********************************************/
/* ENS -- WDM a la Eke Navarro and Steinmetz */
/*********************************************/
double c_ENS(double      M_vir,
	     double      z,
	     cosmo_info *cosmo);
double zbrent_ENS(double func(double,
			      double,
                              double,
			      cosmo_info *,
			      double),
		  double      x1,
		  double      x2,
		  double      tol,
		  double      M_vir,			
                  double      z,
		  cosmo_info *cosmo,
		  double      delviro);
double zbrent_ENS_2(double func(double,
				double,
				cosmo_info *),
		    double      x1,
		    double      x2,
		    double      tol,
		    cosmo_info *cosmo,
		    double      dzc);
double rhoorhoc(double omega,
		double lambda);
double concsuper(double      c,
		 double      M_vir,
                 double      z,
		 cosmo_info *cosmo,
		 double      delviro);
double redcol(double z,double dzc,cosmo_info *cosmo);
double dlnsdlnm(double m,double z,cosmo_info *cosmo);


/*********************************/
/* Correlation function routines */
/*********************************/
double W_k_boxcar(double kR);
double zeta_r(double      R,
	      double      redshift,
	      cosmo_info *cosmo,
	      int         mode,
	      int         component);
double zeta_r_halo(double      M_vir,
		   double      R,
		   double      redshift,
		   cosmo_info *cosmo);
double lk_nonlinear(double      redshift,
                    cosmo_info *cosmo);
double n_effective(double      lk_interp,
                   cosmo_info *cosmo);
double halo_bias_M(double      M_vir,
		   double      z,
		   cosmo_info *cosmo,
		   int         mode,
		   int         component);
double halo_bias_R(double      R,
		   double      z,
		   cosmo_info *cosmo,
		   int         mode,
		   int         component);
double halo_bias(double      M_vir,
		 double      R,
		 double      z,
		 cosmo_info *cosmo,
		 int         mode,
		 int         component);

double Nc_M(double      M,
	    cosmo_info *cosmo);
double Ns_M(double      M,
	    double      z,
	    cosmo_info *cosmo);
double N_M(double      M,
	   double      z,
	   cosmo_info *cosmo);
double halo_density(double      z,
		    cosmo_info *cosmo);
double halo_density_limited(double      lM_max,
			    double      z,
			    cosmo_info *cosmo);
double halo_density_restricted(double       R,
                               double       z,
                               cosmo_info **cosmo);
double zeta_r_halo_1_cs(double      R,
			double      redshift,
			cosmo_info *cosmo);
void init_power_spectrum_halo_1_ss(double       z,
				   cosmo_info **cosmo,
				   int          mode);
double power_spectrum_halo_1_ss(double      k,
				double      z,
				cosmo_info *cosmo,
				int         mode);
double zeta_r_halo_1_ss(double      R,
			double      z,
			cosmo_info *cosmo);
double zeta_r_halo_1(double      R,
		     double      z,
		     cosmo_info *cosmo);
void init_M_lim_restricted(double       z,
			   cosmo_info **cosmo);
double M_lim_restricted(double       R,
			double       z,
			cosmo_info **cosmo);
void init_power_spectrum_halo_2(double       R,
                                double       z,
				cosmo_info **cosmo,
				int          mode);
double power_spectrum_halo_2(double      k,
                             double      R,
			     double      z,
			     cosmo_info *cosmo,
			     int         mode);
double zeta_r_halo_2(double      R,
		     double      z,
		     cosmo_info *cosmo);
#endif
