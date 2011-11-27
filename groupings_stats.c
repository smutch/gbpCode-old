#define  _MAIN
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

int main(int argc, char *argv[]){
  char    filename_in[MAX_FILENAME_LENGTH];
  char    filename_out[MAX_FILENAME_LENGTH];
  char    filename_in_root[MAX_FILENAME_LENGTH];
  int     n_groupings;
  int     i_grouping;
  int     n_halos_per_grouping;
  int     i_halo;
  double *M_halos;
  double *V_halos;
  char   *line=NULL;
  size_t  line_length=0;
  FILE   *fp_in;
  FILE   *fp_out;
  int     i_column=1;
  double  M_min;
  double  M_med;
  double  M_max;
  double  V_min;
  double  V_med;
  double  V_max;

  // Initialization -- MPI etc.
  SID_init(&argc,&argv,NULL);
  strcpy(filename_in_root,    argv[1]);
  n_groupings     =(int) atoi(argv[2]);

  SID_log("Producing statistics for %d groupings...",SID_LOG_OPEN,n_groupings);

  // Write each halo grouping in turn ...
  sprintf(filename_out,"%s_stats.dat",filename_in_root);
  fp_out=fopen(filename_out,"w");
  fprintf(fp_out,"# Stats for %d halo groupings {%s*}\n",n_groupings,filename_in_root);
  fprintf(fp_out,"# Column: (%d) n_halos\n",       i_column++);
  fprintf(fp_out,"#         (%d) M_sub (min)\n",   i_column++);
  fprintf(fp_out,"#         (%d) M_sub (median)\n",i_column++);
  fprintf(fp_out,"#         (%d) M_sub (max)\n",   i_column++);
  fprintf(fp_out,"#         (%d) V_sub (min)\n",   i_column++);
  fprintf(fp_out,"#         (%d) V_sub (median)\n",i_column++);
  fprintf(fp_out,"#         (%d) V_sub (max)\n",   i_column++);
  for(i_grouping=0;i_grouping<n_groupings;i_grouping++){
    SID_log("Analyzing grouping %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_grouping+1,n_groupings);
    sprintf(filename_in,"%s_grouping_%03d.dat",filename_in_root,i_grouping);
    fp_in=fopen(filename_in,"r");
    n_halos_per_grouping=count_lines_data(fp_in);
    M_halos=(double *)SID_malloc(sizeof(double)*n_halos_per_grouping);
    V_halos=(double *)SID_malloc(sizeof(double)*n_halos_per_grouping);
    for(i_halo=0;i_halo<n_halos_per_grouping;i_halo++){
      grab_next_line_data(fp_in,&line,&line_length);
      grab_double(line,1,&(M_halos[i_halo]));
      grab_double(line,2,&(V_halos[i_halo]));
    } 
    fclose(fp_in);

    calc_min(   M_halos,&M_min,n_halos_per_grouping,SID_DOUBLE,CALC_MODE_DEFAULT);
    calc_median(M_halos,&M_med,n_halos_per_grouping,SID_DOUBLE,CALC_MODE_DEFAULT);
    calc_max(   M_halos,&M_max,n_halos_per_grouping,SID_DOUBLE,CALC_MODE_DEFAULT);

    calc_min(   V_halos,&V_min,n_halos_per_grouping,SID_DOUBLE,CALC_MODE_DEFAULT);
    calc_median(V_halos,&V_med,n_halos_per_grouping,SID_DOUBLE,CALC_MODE_DEFAULT);
    calc_max(   V_halos,&V_max,n_halos_per_grouping,SID_DOUBLE,CALC_MODE_DEFAULT);

    fprintf(fp_out,"%d %le %le %le %le %le %le\n",n_halos_per_grouping,M_min,M_med,M_max,V_min,V_med,V_max);

    SID_free(SID_FARG M_halos);
    SID_free(SID_FARG V_halos);
    SID_log("Done.",SID_LOG_CLOSE);
  }
  fclose(fp_out);
  
  // Clean-up
  SID_log("Cleaning-up ...",SID_LOG_OPEN);
  SID_free(SID_FARG line);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
