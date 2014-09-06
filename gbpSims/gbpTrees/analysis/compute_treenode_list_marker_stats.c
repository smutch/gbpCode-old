#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

void compute_treenode_list_marker_stats(tree_info *trees,treenode_list_info *list,tree_markers_info **markers_all,tree_markers_stats_info *stats,int **n_hist_count,int *n_hist){

   int n_stats=3;
   int n_M    =2;
   (*n_hist)  =n_stats+n_M;

   // Allocate histogram arrays -- stats
   int   i_hist=0;
   int **hist  =(int **)SID_calloc(sizeof(int *)*(*n_hist));
   int  *n_bins=(int  *)SID_calloc(sizeof(int)*(*n_hist));
   for(int i_stat=0;i_stat<n_stats;i_stat++){
      n_bins[i_hist]=trees->n_snaps;
      hist[i_hist++]=(int *)SID_calloc(sizeof(int)*n_bins[i_hist]);
   }

   // Allocate histogram arrays -- mass
   double logM_min=7.;
   int    n_M_bins=90;
   double dlogM   =0.1;
   for(int i_M=0;i_M<n_M;i_M++){
      n_bins[i_hist]=n_M_bins;
      hist[i_hist++]=(int *)SID_calloc(sizeof(int)*n_bins[i_hist]);
   }

   // Allocate histogram counts
   if((*n_hist_count)==NULL)
      (*n_hist_count)=(int *)SID_calloc(sizeof(int)*(*n_hist));
   else{
      for(i_hist=0;i_hist<(*n_hist);i_hist++)
         (*n_hist_count)[i_hist]=0;
   }

   // Work with a precompiled markers array if it's been given
   tree_markers_info *markers;
   if(markers_all==NULL)
      markers=(tree_markers_info *)SID_malloc(sizeof(tree_markers_info));

   // Create histograms
   for(int i_list=0;i_list<list->n_list_local;i_list++){
      find_treenode_markers(trees,list->list[i_list],markers_all,markers);
      // Build histogram
      int i_bin;
      for(int i_hist=0;i_hist<(*n_hist);i_hist++){
         i_bin=-1;
         if(i_hist<n_stats){
            tree_node_info *marker=NULL;
            if(i_hist==0)
               marker=markers->half_peak_mass;
            else if(i_hist==1)
               marker=markers->merger_33pc_remnant;
            else if(i_hist==2)
               marker=markers->merger_10pc_remnant;
            else
               SID_trap_error("Behaviour undefined in compute_treenode_list_marker_stats()",ERROR_LOGIC);
            if(marker!=NULL)
               i_bin=marker->snap_tree;
            else
               i_bin=-1;
         }
         else{
            if(i_hist==(n_stats+0)){
               halo_properties_info *properties=fetch_treenode_properties(trees,list->list[i_list]);
               i_bin=(take_log10(properties->M_vir)-logM_min)/dlogM;
            }
            else if(i_hist==(n_stats+1))
               i_bin=(take_log10(markers->M_peak)-logM_min)/dlogM;
            else
               SID_trap_error("Behaviour undefined in compute_treenode_list_marker_stats()",ERROR_LOGIC);
         }
         if(i_bin>=0 && i_bin<n_bins[i_hist]){
            hist[i_hist][i_bin]++;
            (*n_hist_count)[i_hist]++;
         }
      }
   }
   if(markers_all==NULL)
      SID_free(SID_FARG markers);
   for(int i_hist=0;i_hist<(*n_hist);i_hist++){
      SID_Allreduce(SID_IN_PLACE,hist[i_hist],n_bins[i_hist],SID_INT,SID_SUM,SID.COMM_WORLD);
      SID_Allreduce(SID_IN_PLACE,&((*n_hist_count)[i_hist]),1,SID_INT,SID_SUM,SID.COMM_WORLD);
   }

   // Find 68 and 95 percent confidence ranges
   for(int i_hist=0;i_hist<(*n_hist);i_hist++){
      int sum=0;
      for(int i_bin=0;i_bin<n_bins[i_hist];i_bin++)
         sum+=hist[i_hist][i_bin];
      size_t *hist_index=NULL;
      merge_sort(hist[i_hist],(size_t)(n_bins[i_hist]),&hist_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
      int i_peak =hist_index[n_bins[i_hist]-1];
      int i_68_lo=hist_index[n_bins[i_hist]-1];
      int i_68_hi=hist_index[n_bins[i_hist]-1];
      int target =(int)(0.68*(double)sum);
      int accum  =0;
      int i_bin  =0;
      while(accum<=target && i_bin<n_bins[i_hist]){
         size_t idx_i=hist_index[n_bins[i_hist]-i_bin-1];
         if(idx_i<i_68_lo) i_68_lo=idx_i;
         if(idx_i>i_68_hi) i_68_hi=idx_i;
         accum+=hist[i_hist][idx_i];
         i_bin++;
      }
      int i_95_lo=i_68_lo;
      int i_95_hi=i_68_hi;
      target=(int)(0.95*(double)sum);
      while(accum<=target && i_bin<n_bins[i_hist]){
         size_t idx_i=hist_index[n_bins[i_hist]-i_bin-1];
         if(idx_i<i_95_lo) i_95_lo=idx_i;
         if(idx_i>i_95_hi) i_95_hi=idx_i;
         accum+=hist[i_hist][idx_i];
         i_bin++;
      }
      SID_free(SID_FARG hist_index);

      if(i_hist==0){
         stats->t_half_peak_mass_ranges[0][0]=trees->t_list[i_68_lo];
         stats->t_half_peak_mass_ranges[0][1]=trees->t_list[i_68_hi];
         stats->t_half_peak_mass_ranges[1][0]=trees->t_list[i_95_lo];
         stats->t_half_peak_mass_ranges[1][1]=trees->t_list[i_95_hi];
         stats->t_half_peak_mass_peak        =trees->t_list[i_peak];
      }
      else if(i_hist==1){
         stats->t_merger_33pc_ranges[0][0]=trees->t_list[i_68_lo];
         stats->t_merger_33pc_ranges[0][1]=trees->t_list[i_68_hi];
         stats->t_merger_33pc_ranges[1][0]=trees->t_list[i_95_lo];
         stats->t_merger_33pc_ranges[1][1]=trees->t_list[i_95_hi];
         stats->t_merger_33pc_peak        =trees->t_list[i_peak];
      }
      else if(i_hist==2){
         stats->t_merger_10pc_ranges[0][0]=trees->t_list[i_68_lo];
         stats->t_merger_10pc_ranges[0][1]=trees->t_list[i_68_hi];
         stats->t_merger_10pc_ranges[1][0]=trees->t_list[i_95_lo];
         stats->t_merger_10pc_ranges[1][1]=trees->t_list[i_95_hi];
         stats->t_merger_10pc_peak        =trees->t_list[i_peak];
      }
      else if(i_hist==3){
         stats->M_peak_ranges[0][0]=logM_min+(double)(i_68_lo)*dlogM;
         stats->M_peak_ranges[0][1]=logM_min+(double)(i_68_hi)*dlogM;
         stats->M_peak_ranges[1][0]=logM_min+(double)(i_95_lo)*dlogM;
         stats->M_peak_ranges[1][1]=logM_min+(double)(i_95_hi)*dlogM;
         stats->M_peak_peak        =logM_min+(double)( i_peak)*dlogM;
      }
      else if(i_hist==4){
         stats->M_vir_ranges[0][0]=logM_min+(double)(i_68_lo)*dlogM;
         stats->M_vir_ranges[0][1]=logM_min+(double)(i_68_hi)*dlogM;
         stats->M_vir_ranges[1][0]=logM_min+(double)(i_95_lo)*dlogM;
         stats->M_vir_ranges[1][1]=logM_min+(double)(i_95_hi)*dlogM;
         stats->M_vir_peak        =logM_min+(double)( i_peak)*dlogM;
      }
      else
         SID_trap_error("Behaviour undefined in compute_treenode_list_marker_stats()",ERROR_LOGIC);
   }

   // Clean-up
   for(int i_hist=0;i_hist<(*n_hist);i_hist++)
      SID_free(SID_FARG hist[i_hist]);
   SID_free(SID_FARG hist);
   SID_free(SID_FARG n_bins);

}

