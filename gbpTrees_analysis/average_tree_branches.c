#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void average_tree_branches(const char *catalog_name){
  SID_log("Processing tree tracks in catalog {%s}...",SID_LOG_OPEN,catalog_name);

  // Master Rank does all the work
  FILE *fp_tracks_in=NULL;
  if(SID.I_am_Master){
     // Create and open the output files
     char   filename_tracks_out[MAX_FILENAME_LENGTH];
     sprintf(filename_tracks_out,"%s_tracks.dat",catalog_name);
     fp_tracks_in=fopen(filename_tracks_out,"r");

     // Write header for tracks file
     int n_list;
     int n_snaps;
     fread(&n_list, sizeof(int),1,fp_tracks_in);
     fread(&n_snaps,sizeof(int),1,fp_tracks_in);
     int    *snap_list=(int    *)SID_malloc(sizeof(int)   *n_snaps);
     double *z_list   =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *t_list   =(double *)SID_malloc(sizeof(double)*n_snaps);
     fread(snap_list,sizeof(int),   n_snaps,fp_tracks_in);
     fread(z_list,   sizeof(double),n_snaps,fp_tracks_in);
     fread(t_list,   sizeof(double),n_snaps,fp_tracks_in);

     // Allocate some temporary arrays for the tracks
     double   M_min   = 6.;
     double   M_max   =16.;
     int      n_M_bins=200;
     double   dM      =(M_max-M_min)/(double)n_M_bins;
     double   inv_dM  =1./dM;
     int    **M_hist  =(int **)SID_malloc(sizeof(int *)*n_snaps);
     for(int i_snap=0;i_snap<n_snaps;i_snap++)
        M_hist[i_snap]=(int *)SID_calloc(sizeof(int)*n_M_bins);
     int    *i_z_track=(int    *)SID_malloc(sizeof(int)*n_snaps);;
     int    *idx_track=(int    *)SID_malloc(sizeof(int)*n_snaps);;
     double *M_track  =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *x_track  =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *y_track  =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *z_track  =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *vx_track =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *vy_track =(double *)SID_malloc(sizeof(double)*n_snaps);
     double *vz_track =(double *)SID_malloc(sizeof(double)*n_snaps);

     // Process each track in turn
     for(int i_list=0;i_list<n_list;i_list++){
        int n_track;
        // Read track
        fread(&n_track, sizeof(int),   1,      fp_tracks_in);
        fread(i_z_track,sizeof(int),   n_track,fp_tracks_in);
        fread(idx_track,sizeof(int),   n_track,fp_tracks_in);
        fread(x_track,  sizeof(double),n_track,fp_tracks_in);
        fread(y_track,  sizeof(double),n_track,fp_tracks_in);
        fread(z_track,  sizeof(double),n_track,fp_tracks_in);
        fread(vx_track, sizeof(double),n_track,fp_tracks_in);
        fread(vy_track, sizeof(double),n_track,fp_tracks_in);
        fread(vz_track, sizeof(double),n_track,fp_tracks_in);
        fread(M_track,  sizeof(double),n_track,fp_tracks_in);
        // Build the M-histograms
        for(int i_track=0;i_track<n_track;i_track++){
           int i_bin=(int)((take_log10(M_track[i_track])-M_min)*inv_dM);
           if(i_bin>=0 && i_bin<n_M_bins)
              M_hist[i_z_track[i_track]][i_bin]++;
        }
     } // for i_list
     fclose(fp_tracks_in);

     // Build confidence intervals for M-track
     int    *n_i    =(int    *)SID_calloc(sizeof(int)*n_snaps);
     double *M_peak =(double *)SID_calloc(sizeof(double)*n_snaps);
     double *M_68_lo=(double *)SID_calloc(sizeof(double)*n_snaps);
     double *M_68_hi=(double *)SID_calloc(sizeof(double)*n_snaps);
     for(int i_snap=0;i_snap<n_snaps;i_snap++){
        size_t *M_hist_index=NULL;
        for(int i_bin=0;i_bin<n_M_bins;i_bin++)
           n_i[i_snap]+=M_hist[i_snap][i_bin];
        if(n_i[i_snap]>0){
           merge_sort(M_hist[i_snap],(size_t)n_M_bins,&M_hist_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
           int i_peak =M_hist_index[n_M_bins-1];
           int i_68_lo=M_hist_index[n_M_bins-1];
           int i_68_hi=M_hist_index[n_M_bins-1];
           int target =(int)(0.68*(double)n_i[i_snap]);
           int accum  =0;
           int i_bin  =0;
           while(accum<=target && i_bin<n_M_bins){
              size_t idx_i=M_hist_index[n_M_bins-i_bin-1];
              if(idx_i<i_68_lo) i_68_lo=idx_i;
              if(idx_i>i_68_hi) i_68_hi=idx_i;
              accum+=M_hist[i_snap][idx_i];
              i_bin++;
           }
           M_peak[i_snap] =M_min+((double)i_peak +0.5)*dM;
           M_68_lo[i_snap]=M_min+((double)i_68_lo+0.5)*dM;
           M_68_hi[i_snap]=M_min+((double)i_68_hi+0.5)*dM;
           SID_free(SID_FARG M_hist_index);
        }
        else{
           M_peak[i_snap] =-1;
           M_68_lo[i_snap]=-1;
           M_68_hi[i_snap]=-1;
        }
     }

     // Write results
     char  filename_out[MAX_FILENAME_LENGTH];
     FILE *fp_out;
     sprintf(filename_out,"%s_tracks.ascii",catalog_name);
     fp_out=fopen(filename_out,"w");
     for(int i_snap=0;i_snap<n_snaps;i_snap++)
        fprintf(fp_out,"%le %le %d %le %le %le\n",
                       z_list[i_snap],
                       t_list[i_snap]/S_PER_YEAR,
                       n_i[i_snap],
                       M_peak[i_snap],
                       M_68_lo[i_snap],
                       M_68_hi[i_snap]);
     fclose(fp_out);

     // Clean-up
     for(int i_snap=0;i_snap<n_snaps;i_snap++)
        SID_free(SID_FARG M_hist[i_snap]);
     SID_free(SID_FARG M_hist);
     SID_free(SID_FARG n_i);
     SID_free(SID_FARG M_peak);
     SID_free(SID_FARG M_68_lo);
     SID_free(SID_FARG M_68_hi);
     SID_free(SID_FARG i_z_track);
     SID_free(SID_FARG idx_track);
     SID_free(SID_FARG M_track);
     SID_free(SID_FARG x_track);
     SID_free(SID_FARG y_track);
     SID_free(SID_FARG z_track);
     SID_free(SID_FARG vx_track);
     SID_free(SID_FARG vy_track);
     SID_free(SID_FARG vz_track);
     SID_free(SID_FARG snap_list);
     SID_free(SID_FARG z_list);
     SID_free(SID_FARG t_list);
  } // if I_am_Master
  SID_Barrier(SID.COMM_WORLD);

  SID_log("Done.",SID_LOG_CLOSE);
}

