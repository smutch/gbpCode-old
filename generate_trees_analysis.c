#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void bin_progenitor_ages_recursive_local(tree_info       *trees,
                                         tree_node_info  *descendant,
                                         double          *t_start_mpr,
                                         int             *n_particles_range,
                                         int              n_ranges,
                                         double          *hist_lo,
                                         double          *hist_hi,
                                         int            **hist_all,
                                         int            **hist_mpr,
                                         int              n_hist){
   tree_node_info *main_progenitor;
   double          mpr_t_start;
   double          t_descendant;
   t_descendant   =trees->t_list[descendant->snap_tree]/(1e9*S_PER_YEAR);
   main_progenitor=descendant->progenitor_first;
   if(main_progenitor!=NULL){
      int    n_halos=1;
      
      // Process the main progenitor first
      bin_progenitor_ages_recursive_local(trees,
                                          main_progenitor,
                                          &mpr_t_start,
                                          n_particles_range,
                                          n_ranges,
                                          hist_lo,
                                          hist_hi,
                                          hist_all,
                                          hist_mpr,
                                          n_hist);
      if(descendant->descendant==NULL){
         double delta_t=t_descendant-mpr_t_start;
         int    flag_continue=TRUE;
         int    i_hist;
         int    i_range;
         for(i_hist=0;i_hist<n_hist && flag_continue;i_hist++){
            if(delta_t>=hist_lo[i_hist] && delta_t<hist_hi[i_hist]){
               for(i_range=0;i_range<n_ranges;i_range++){
                  if(main_progenitor->n_particles>n_particles_range[i_range]){
                     hist_mpr[i_range][i_hist]++;
                     hist_all[i_range][i_hist]++;
                  }
               }
               flag_continue=FALSE;
            }
         }
      }

      // Then process the others
      tree_node_info *current_progenitor;
      current_progenitor=main_progenitor->progenitor_next;
      while(current_progenitor!=NULL){
         double cpr_t_start;
         if(n_halos>=descendant->n_progenitors)
            SID_trap_error("Progenitor count exceeded in compute_progenitor_order_recursive().",ERROR_LOGIC);
         bin_progenitor_ages_recursive_local(trees,
                                             current_progenitor,
                                             &cpr_t_start,
                                             n_particles_range,
                                             n_ranges,
                                             hist_lo,
                                             hist_hi,
                                             hist_all,
                                             hist_mpr,
                                             n_hist);
         double delta_t=t_descendant-cpr_t_start;
         int    flag_continue=TRUE;
         int    i_hist;
         int    i_range;
         int    flag_fragmented=check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_FRAGMENTED_RETURNED)||
                                check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_FRAGMENTED_EXCHANGED);
         if(!flag_fragmented){
            for(i_hist=0;i_hist<n_hist && flag_continue;i_hist++){
               if(delta_t>=hist_lo[i_hist] && delta_t<hist_hi[i_hist]){
                  for(i_range=0;i_range<n_ranges;i_range++){
                     if(current_progenitor->n_particles>n_particles_range[i_range])
                        hist_all[i_range][i_hist]++;
                  }
                  flag_continue=FALSE;
               }
            }
         }

         n_halos++;
         current_progenitor=current_progenitor->progenitor_next;
      }
   }

   if(t_start_mpr!=NULL){
      if(descendant->progenitor_first==NULL)
         (*t_start_mpr)=t_descendant;
      else
         (*t_start_mpr)=mpr_t_start;
   }

}

void generate_trees_analysis(tree_info *trees,char *filename_out_root){

  // Compute merger rates ...
  SID_log("Performing analysis of trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  int   i_snap;
  int   i_range;
  int   n_ranges;

  // ... set particle ranges ...
  int *n_particles_range;
  n_particles_range           =(int *)SID_malloc(sizeof(int)*n_ranges);
  i_range                     =0;
  n_ranges                    =5;
  n_particles_range[i_range++]=0;
  n_particles_range[i_range++]=64;
  for(;i_range<n_ranges;i_range++)
     n_particles_range[i_range]=n_particles_range[i_range-1]*8;

  // ... create merger and halo count arrays ...
  int **n_mergers;
  int **n_halos;
  n_mergers=(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_halos  =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     n_mergers[i_snap]=(int *)SID_calloc(sizeof(int)*n_ranges);
     n_halos[i_snap]  =(int *)SID_calloc(sizeof(int)*n_ranges);
  }

  // ... loop over each snapshot ...
  SID_log("Counting mergers...",SID_LOG_OPEN);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;

     // ... loop over each halo for this snapshot ...
     current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        tree_node_info *main_progenitor;
        tree_node_info *merging_progenitor;

        // ... get the main progenitor ...
        main_progenitor   =current_halo->progenitor_first;

        // ... loop over each merging progenitor for this halo ...
        if(main_progenitor!=NULL){
           merging_progenitor=main_progenitor->progenitor_next;
           while(merging_progenitor!=NULL){
              // ... process each merger ...
              int    flag_fragmented=check_mode_for_flag(merging_progenitor->tree_case,TREE_CASE_FRAGMENTED_RETURNED)||
                                     check_mode_for_flag(merging_progenitor->tree_case,TREE_CASE_FRAGMENTED_EXCHANGED);
              if(!flag_fragmented){
                 for(i_range=0;i_range<n_ranges;i_range++){
                    if(merging_progenitor->n_particles>n_particles_range[i_range])
                       n_mergers[i_snap][i_range]++;
                 }
              }
              merging_progenitor=merging_progenitor->progenitor_next;
           }
        }

        // ... count the number of halos ...
        for(i_range=0;i_range<n_ranges;i_range++){
           if(current_halo->n_particles>n_particles_range[i_range])
              n_halos[i_snap][i_range]++;
        }

        current_halo=current_halo->next_neighbour;
     }
     SID_Allreduce(SID_IN_PLACE,n_mergers[i_snap],n_ranges,SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,n_halos[i_snap],  n_ranges,SID_INT,SID_SUM,SID.COMM_WORLD);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write mergers file
  char  filename_out_mergers[MAX_FILENAME_LENGTH];
  FILE *fp;
  int   i_column=1;
  sprintf(filename_out_mergers,"%s_mergers.txt",filename_out_root);
  SID_log("Writing merger statistics to {%s}...",SID_LOG_OPEN,filename_out_mergers);
  fp=fopen(filename_out_mergers,"w");
  fprintf(fp,"# Merger statistics for {%s}\n",filename_out_root);
  fprintf(fp,"# Column (%02d): Snapshot number\n",     i_column++);
  fprintf(fp,"#        (%02d): Snapshot redshift\n",   i_column++);
  fprintf(fp,"#        (%02d): Snapshot time [Gyrs]\n",i_column++);
  for(i_range=0;i_range<n_ranges;i_range++){
     fprintf(fp,"#        (%02d): No. of halos;              >=%04d particles\n",i_column++,n_particles_range[i_range]);
     fprintf(fp,"#        (%02d): No. of mergers;            >=%04d particles\n",i_column++,n_particles_range[i_range]);
     fprintf(fp,"#        (%02d): Fraction of halos merging; >=%04d particles\n",i_column++,n_particles_range[i_range]);
  }
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     fprintf(fp,"%3d %10.3lf %10.3lf",
             trees->snap_list[i_snap],
             trees->z_list[i_snap],
             trees->t_list[i_snap]/(1e9*S_PER_YEAR));
     for(i_range=0;i_range<n_ranges;i_range++){
        double fraction;
        if(n_halos[i_snap][i_range]>0)
           fraction=(double)n_mergers[i_snap][i_range]/(double)n_halos[i_snap][i_range];
        else
           fraction=0.;
        fprintf(fp," %5d %5d %10.4lf",
                   n_halos[i_snap][i_range],
                   n_mergers[i_snap][i_range],
                   fraction);
     }
     fprintf(fp,"\n");
  }
  fclose(fp);
  SID_log("Done.",SID_LOG_CLOSE);

  // Creating a histogram of progenitor ages
  int    i_hist;
  int    n_hist      = 26;
  double hist_min    = 0.;
  double hist_max    =13.;
  double (*hist_lo)  =(double  *)SID_malloc(sizeof(double)*n_hist);
  double (*hist_hi)  =(double  *)SID_malloc(sizeof(double)*n_hist);
  int    (**hist_all)=(int    **)SID_malloc(sizeof(int *)*n_ranges);
  int    (**hist_mpr)=(int    **)SID_malloc(sizeof(int *)*n_ranges);
  double hist_step =(hist_max-hist_min)/(double)n_hist;
  for(i_range=0;i_range<n_ranges;i_range++){
    hist_all[i_range]=(int *)SID_calloc(sizeof(int)*n_hist);
    hist_mpr[i_range]=(int *)SID_calloc(sizeof(int)*n_hist);
  }
  for(i_hist=0;i_hist<n_hist;i_hist++){
     if(i_hist==0)
        hist_lo[i_hist]=0.;
     else
        hist_lo[i_hist]=hist_hi[i_hist-1];
     hist_hi[i_hist]=hist_lo[i_hist]+hist_step;
  }
  hist_hi[n_hist-1]=hist_max;
     
  SID_log("Counting progenitor ages...",SID_LOG_OPEN|SID_LOG_TIMER);
  // ... loop over each halo for this snapshot ...
  tree_node_info *current_halo=trees->first_neighbour_subgroups[trees->n_snaps-1];
  while(current_halo!=NULL){
     bin_progenitor_ages_recursive_local(trees,
                                         current_halo,
                                         NULL,
                                         n_particles_range,
                                         n_ranges,
                                         hist_lo,
                                         hist_hi,
                                         hist_all,
                                         hist_mpr,
                                         n_hist);
     current_halo=current_halo->next_neighbour;
  }
  for(i_range=0;i_range<n_ranges;i_range++){
     SID_Allreduce(SID_IN_PLACE,hist_all[i_range],n_hist,SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,hist_mpr[i_range],n_hist,SID_INT,SID_SUM,SID.COMM_WORLD);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write progenitor ages
  char  filename_out_ages[MAX_FILENAME_LENGTH];
  i_column=1;
  sprintf(filename_out_ages,"%s_ages.txt",filename_out_root);
  SID_log("Writing progenitor ages to {%s}...",SID_LOG_OPEN,filename_out_ages);
  fp=fopen(filename_out_ages,"w");
  fprintf(fp,"# Merger statistics for {%s}\n",filename_out_root);
  fprintf(fp,"# Column (%02d): Age bin (lower limit) [Gyrs]\n", i_column++);
  fprintf(fp,"# Column (%02d): Age bin (upper limit) [Gyrs]\n", i_column++);
  for(i_range=0;i_range<n_ranges;i_range++){
     fprintf(fp,"#        (%02d): No. of main progenitors;  >=%04d particles\n",i_column++,n_particles_range[i_range]);
     fprintf(fp,"#        (%02d): No. of progenitors (all); >=%04d particles\n",i_column++,n_particles_range[i_range]);
  }
  for(i_hist=0;i_hist<n_hist;i_hist++){
     fprintf(fp,"%10.3lf %10.3lf",hist_lo[i_hist],hist_hi[i_hist]);
     for(i_range=0;i_range<n_ranges;i_range++)
        fprintf(fp," %10.4lf %10.4lf",
                   (double)hist_mpr[i_range][i_hist],
                   (double)hist_all[i_range][i_hist]);
     fprintf(fp,"\n");
  }
  fclose(fp);
  for(i_range=0;i_range<n_ranges;i_range++){
     SID_free(SID_FARG hist_mpr[i_range]);
     SID_free(SID_FARG hist_all[i_range]);
  }
  SID_free(SID_FARG hist_mpr);
  SID_free(SID_FARG hist_all);
  SID_free(SID_FARG hist_lo);
  SID_free(SID_FARG hist_hi);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     SID_free(SID_FARG n_mergers[i_snap]);
     SID_free(SID_FARG n_halos[i_snap]);
  }
  SID_free(SID_FARG n_mergers);
  SID_free(SID_FARG n_halos);
  SID_free(SID_FARG n_particles_range);

  SID_log("Done.",SID_LOG_CLOSE);
}

