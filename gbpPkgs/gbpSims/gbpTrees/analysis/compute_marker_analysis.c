#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpHighZ.h>

void compute_marker_analysis(tree_info *trees,tree_markers_info **markers,hist_info *M_hist,hist_info *xoff_hist,hist_info *SSFctn_hist,const char *catalog_root,const char *filename_out_root,int mode){

  // Compute merger rates ...
  char filename_out_root_group[MAX_FILENAME_LENGTH];
  if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_GROUPS)){
     SID_log("Performing GROUP marker analysis...",SID_LOG_OPEN|SID_LOG_TIMER);
     sprintf(filename_out_root_group,"%s_groups",filename_out_root);
  }
  else if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS)){
     SID_log("Performing SUBGROUP marker analysis...",SID_LOG_OPEN|SID_LOG_TIMER);
     sprintf(filename_out_root_group,"%s_subgroups",filename_out_root);
  }
  else
     SID_trap_error("group/subgroup mode has not been properly specified in compute_marker_analysis().",ERROR_LOGIC);


  char  filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"%s.txt",filename_out_root_group);
  treenode_list_info **lists=(treenode_list_info **)SID_malloc(sizeof(treenode_list_info *)*n_bins);

  int i_bin;
  for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
     SID_log("Processing snapshot #%03d of %03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_snap+1,trees->n_snaps);

     // Count halos
     tree_node_info *current_halo=NULL;
     if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_GROUPS))
        current_halo=trees->first_neighbour_groups[i_snap];
     else if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS))
        current_halo=trees->first_neighbour_subgroups[i_snap];
     for(i_bin=0;i_bin<n_bins;i_bin++)
        n_bin[i_bin]=0;
     while(current_halo!=NULL){
        double logM_i=take_log10(fetch_treenode_M_vir(trees,current_halo));
        int    bin   =(int)((logM_i-logM_min)/dlogM);
        if(bin>=0 && bin<n_bins)
           n_bin[bin]++;
        current_halo=current_halo->next_neighbour;
     }

     // Create lists of halos for each mass bin
     for(i_bin=0;i_bin<n_bins;i_bin++){
        char list_name[32];
        sprintf(list_name,"bin_%03d",i_bin);
        init_treenode_list(list_name,n_bin[i_bin],&(lists[i_bin]));
     }
     if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_GROUPS))
        current_halo=trees->first_neighbour_groups[i_snap];
     else if(check_mode_for_flag(mode,COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS))
        current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        double logM_i=take_log10(fetch_treenode_M_vir(trees,current_halo));
        int    bin   =(int)((logM_i-logM_min)/dlogM);
        if(bin>=0 && bin<n_bins)
           add_to_treenode_list(lists[bin],current_halo);
        current_halo=current_halo->next_neighbour;
     }

     // Finalize lists
     for(i_bin=0;i_bin<n_bins;i_bin++)
        finalize_treenode_list(trees,lists[i_bin]);

     // Write header
     FILE *fp_out=NULL;
     if(SID.I_am_Master){
        if(i_snap==0){
           fp_out=fopen(filename_out,"w");
           int i_column=1;
           fprintf(fp_out,"#  Mass accretion statistics for {%s}\n",catalog_root);
           fprintf(fp_out,"#  Column (%03d): Snapshot\n",    i_column++);
           fprintf(fp_out,"#         (%03d): Redshift\n",    i_column++);
           fprintf(fp_out,"#         (%03d): t_snap [yrs]\n",i_column++);
           for(i_bin=0;i_bin<n_bins;i_bin++){
              double logM_lo=logM_min+dlogM*(double)(i_bin);
              double logM_hi=logM_min+dlogM*(double)(i_bin+1);
              fprintf(fp_out,"#         (%03d): n_halos [log_M=%5.2lf to %5.2lf]\n",i_column++,logM_lo,logM_hi);
              fprintf(fp_out,"#         (%03d): n_t(half_peak_mass_peak)\n",                         i_column++);
              fprintf(fp_out,"#         (%03d): t(half_peak_mass_peak) [yrs]\n",                     i_column++);
              fprintf(fp_out,"#         (%03d): t(half_peak_mass_peak) [yrs] - 68%% confidence lo\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(half_peak_mass_peak) [yrs] - 68%% confidence hi\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(half_peak_mass_peak) [yrs] - 95%% confidence lo\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(half_peak_mass_peak) [yrs] - 95%% confidence hi\n",i_column++);
              fprintf(fp_out,"#         (%03d): n_t(last 3:1 merger)\n",                             i_column++);
              fprintf(fp_out,"#         (%03d): t(last  3:1 merger)    [yrs]\n",                     i_column++);
              fprintf(fp_out,"#         (%03d): t(last  3:1 merger)    [yrs] - 68%% confidence lo\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(last  3:1 merger)    [yrs] - 68%% confidence hi\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(last  3:1 merger)    [yrs] - 95%% confidence lo\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(last  3:1 merger)    [yrs] - 95%% confidence hi\n",i_column++);
              fprintf(fp_out,"#         (%03d): n_t(last 10:1 merger)\n",                            i_column++);
              fprintf(fp_out,"#         (%03d): t(last 10:1 merger)    [yrs]\n",                     i_column++);
              fprintf(fp_out,"#         (%03d): t(last 10:1 merger)    [yrs] - 68%% confidence lo\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(last 10:1 merger)    [yrs] - 68%% confidence hi\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(last 10:1 merger)    [yrs] - 95%% confidence lo\n",i_column++);
              fprintf(fp_out,"#         (%03d): t(last 10:1 merger)    [yrs] - 95%% confidence hi\n",i_column++);
              fprintf(fp_out,"#         (%03d): n_M_vir\n",                                          i_column++);
              fprintf(fp_out,"#         (%03d): M_vir(z) [h^{-1} M_sol]\n",                          i_column++);
              fprintf(fp_out,"#         (%03d): M_vir(z) [h^{-1} M_sol] - 68%% confidence lo\n",     i_column++);
              fprintf(fp_out,"#         (%03d): M_vir(z) [h^{-1} M_sol] - 68%% confidence hi\n",     i_column++);
              fprintf(fp_out,"#         (%03d): M_vir(z) [h^{-1} M_sol] - 95%% confidence lo\n",     i_column++);
              fprintf(fp_out,"#         (%03d): M_vir(z) [h^{-1} M_sol] - 95%% confidence hi\n",     i_column++);
              fprintf(fp_out,"#         (%03d): n_M_peak\n",                                         i_column++);
              fprintf(fp_out,"#         (%03d): M_peak(z) [h^{-1} M_sol]\n",                         i_column++);
              fprintf(fp_out,"#         (%03d): M_peak(z) [h^{-1} M_sol] - 68%% confidence lo\n",    i_column++);
              fprintf(fp_out,"#         (%03d): M_peak(z) [h^{-1} M_sol] - 68%% confidence hi\n",    i_column++);
              fprintf(fp_out,"#         (%03d): M_peak(z) [h^{-1} M_sol] - 95%% confidence lo\n",    i_column++);
              fprintf(fp_out,"#         (%03d): M_peak(z) [h^{-1} M_sol] - 95%% confidence hi\n",    i_column++);
           }
        }
        else
           fp_out=fopen(filename_out,"a");
        fprintf(fp_out,"%3d %le %le",trees->snap_list[i_snap],trees->z_list[i_snap],trees->t_list[i_snap]/S_PER_YEAR);
     }

     // Compute and write statistics
     int *n_hist_count=NULL;
     for(i_bin=0;i_bin<n_bins;i_bin++){
        tree_markers_stats_info stats;
        int                     n_hist;
        compute_treenode_list_marker_stats(trees,lists[i_bin],markers,&stats,&n_hist_count,&n_hist);
        if(SID.I_am_Master){
           fprintf(fp_out," %d %d %le %le %le %le %le %d %le %le %le %le %le %d %le %le %le %le %le %d %le %le %le %le %le %d %le %le %le %le %le",
                          lists[i_bin]->n_list,
                          n_hist_count[0],
                          stats.t_half_peak_mass_peak/S_PER_YEAR,
                          stats.t_half_peak_mass_ranges[0][0]/S_PER_YEAR,
                          stats.t_half_peak_mass_ranges[0][1]/S_PER_YEAR,
                          stats.t_half_peak_mass_ranges[1][0]/S_PER_YEAR,
                          stats.t_half_peak_mass_ranges[1][1]/S_PER_YEAR,
                          n_hist_count[1],
                          stats.t_merger_33pc_peak/S_PER_YEAR,
                          stats.t_merger_33pc_ranges[0][0]/S_PER_YEAR,
                          stats.t_merger_33pc_ranges[0][1]/S_PER_YEAR,
                          stats.t_merger_33pc_ranges[1][0]/S_PER_YEAR,
                          stats.t_merger_33pc_ranges[1][1]/S_PER_YEAR,
                          n_hist_count[2],
                          stats.t_merger_10pc_peak/S_PER_YEAR,
                          stats.t_merger_10pc_ranges[0][0]/S_PER_YEAR,
                          stats.t_merger_10pc_ranges[0][1]/S_PER_YEAR,
                          stats.t_merger_10pc_ranges[1][0]/S_PER_YEAR,
                          stats.t_merger_10pc_ranges[1][1]/S_PER_YEAR,
                          n_hist_count[3],
                          stats.M_vir_peak,
                          stats.M_vir_ranges[0][0],
                          stats.M_vir_ranges[0][1],
                          stats.M_vir_ranges[1][0],
                          stats.M_vir_ranges[1][1],
                          n_hist_count[4],
                          stats.M_peak_peak,
                          stats.M_peak_ranges[0][0],
                          stats.M_peak_ranges[0][1],
                          stats.M_peak_ranges[1][0],
                          stats.M_peak_ranges[1][1]);
        }
     }
     SID_free(SID_FARG n_hist_count);
     if(SID.I_am_Master){
        fprintf(fp_out,"\n");
        fclose(fp_out);
     }

     // Clean-up
     for(i_bin=0;i_bin<n_bins;i_bin++)
        free_treenode_list(&(lists[i_bin]));

     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up
  SID_free(SID_FARG n_bin);
  SID_free(SID_FARG lists);

  SID_log("Done.",SID_LOG_CLOSE);
}

