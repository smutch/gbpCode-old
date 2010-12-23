#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpRNG.h>
#include <gbpMCMC.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>

void compute_MCMC_chain_stats(double **x,int n_x,int n_avg,double *x_min,double *x_bar_in,double *x_max,double *x_sigma,double **auto_cor,double *slopes,double *dP_sub,double *Chi2_in,double *Chi2_min,double *Chi2_avg,double *Chi2_max){
  int     i_P;
  int     i_avg;
  double *x_bar;
  double *s_temp;
  double  c0,c1,cov00,cov01,cov11,sumsq;
  double  covar;
  double  var;

  // Allocate a temporary array for the mean if it hasn't been given/isn't requested but is needed
  if((x_sigma!=NULL || auto_cor!=NULL || slopes!=NULL) && x_bar_in==NULL)
    x_bar=(double *)SID_malloc(sizeof(double)*n_x);
  else
    x_bar=x_bar_in;

  // Compute the mean (if it's requested/needed)
  if(x_bar!=NULL){
    for(i_P=0;i_P<n_x;i_P++){
      x_bar[i_P]=0.;
      for(i_avg=0;i_avg<n_avg;i_avg++)
        x_bar[i_P]+=x[i_P][i_avg];
      x_bar[i_P]/=(double)n_avg;
    }
  }
  
  // Compute the std deviation (if it's requested/needed)
  if(x_sigma!=NULL){
    for(i_P=0;i_P<n_x;i_P++){
      x_sigma[i_P]=0.;
      for(i_avg=0;i_avg<n_avg;i_avg++)
        x_sigma[i_P]+=pow(x[i_P][i_avg]-x_bar[i_P],2.);
      x_sigma[i_P]=sqrt(x_sigma[i_P]/(double)n_avg);
    }
  }

  // Compute the slope across the interval (if it's requested/needed)
  if(slopes!=NULL){
    s_temp=(double *)SID_malloc(sizeof(double)*n_avg);
    for(i_avg=0;i_avg<n_avg;i_avg++)
      s_temp[i_avg]=(double)i_avg;
    for(i_P=0;i_P<n_x;i_P++){
      gsl_fit_linear(s_temp,1,x[i_P],1,n_avg,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
      slopes[i_P]=c1;
      if(dP_sub!=NULL)
        dP_sub[i_P]=sqrt(sumsq/(double)n_avg);
    }
    SID_free(SID_FARG s_temp);
  }
  
  // Compute the auto correlation function (if it's requested/needed)
  if(auto_cor!=NULL){
    for(i_P=0;i_P<n_x;i_P++){
      for(i_P=1;i_P<n_avg;i_P++){
        for(i_avg=0,var=0.,covar=0.;i_avg<(n_avg-i_P);i_avg++){
          covar+=(x[i_P][i_avg]-x_bar[i_P])*(x[i_P][i_avg+i_P]-x_bar[i_P]);
          var  +=(x[i_P][i_avg]-x_bar[i_P])*(x[i_P][  i_avg  ]-x_bar[i_P]);
        }
        auto_cor[i_P][i_P-1]=covar/var;
      }
    }
  }

  // Compute Chi2 statistics
  if(Chi2_in!=NULL){
    (*Chi2_min)=Chi2_in[0];
    (*Chi2_avg)=Chi2_in[0];
    (*Chi2_max)=Chi2_in[0];
    for(i_avg=1;i_avg<n_avg;i_avg++){
      (*Chi2_min) =MIN((*Chi2_min),Chi2_in[i_avg]);
      (*Chi2_avg)+=Chi2_in[i_avg];
      (*Chi2_max) =MAX((*Chi2_max),Chi2_in[i_avg]);
    }
    (*Chi2_avg)/=(double)n_avg;
  }

  if(x_bar!=x_bar_in)
    SID_free(SID_FARG x_bar);
}

