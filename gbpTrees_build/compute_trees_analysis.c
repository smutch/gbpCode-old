#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <assert.h>

void compute_trees_analysis(tree_info *trees,char *filename_out_root){

  float **match_score_local=(float **)ADaPS_fetch(trees->data,"match_score_subgroups");

  // Compute merger rates ...
  SID_log("Performing analysis of trees...",SID_LOG_OPEN|SID_LOG_TIMER);
  int   i_snap;
  int   i_np;
  int   i_xi;

  // Set output directory
  char filename_out_root_dir[MAX_FILENAME_LENGTH];
  char filename_out_dir[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root_dir,"%s_tree_analysis/",filename_out_root);
  mkdir(filename_out_root_dir,02755);

  // Create halo size bins
  int  i_bin     =0;
  int  n_np      =18;
  int *bin_np_lo =(int *)SID_malloc(sizeof(int)*n_np);
  int *bin_np_hi =(int *)SID_malloc(sizeof(int)*n_np);
  bin_np_lo[i_bin]=32;
  bin_np_hi[i_bin]=bin_np_lo[i_bin]*2;
  for(i_bin++;i_bin<n_np;i_bin++){
     bin_np_lo[i_bin]=bin_np_hi[i_bin-1];
     bin_np_hi[i_bin]=bin_np_lo[i_bin]*2;
  }

  // Create merger and halo count arrays
  int   n_bins=n_np*n_np;
  int **n_halos;
  int **n_dropped;
  int **n_emerged;
  int **n_fragmented;
  int **n_mergers;
  int **n_emerged_ages;
  int **n_dropped_ages;
  int **n_fragmented_ages;
  int **n_mergers_ages;
  n_halos          =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_emerged        =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_dropped        =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_fragmented     =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_mergers        =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_emerged_ages   =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_dropped_ages   =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_fragmented_ages=(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  n_mergers_ages   =(int **)SID_malloc(sizeof(int *)*trees->n_snaps);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     n_halos[i_snap]          =(int *)SID_calloc(sizeof(int)*n_np);
     n_emerged[i_snap]        =(int *)SID_calloc(sizeof(int)*n_np);
     n_dropped[i_snap]        =(int *)SID_calloc(sizeof(int)*n_np);
     n_fragmented[i_snap]     =(int *)SID_calloc(sizeof(int)*n_bins);
     n_mergers[i_snap]        =(int *)SID_calloc(sizeof(int)*n_bins);
     n_emerged_ages[i_snap]   =(int *)SID_calloc(sizeof(int)*n_np*i_snap);
     n_dropped_ages[i_snap]   =(int *)SID_calloc(sizeof(int)*n_np*i_snap);
     n_fragmented_ages[i_snap]=(int *)SID_calloc(sizeof(int)*n_bins*i_snap);
     n_mergers_ages[i_snap]   =(int *)SID_calloc(sizeof(int)*n_bins*i_snap);
  }

  // Loop over each snapshot ...
  SID_log("Performing counts...",SID_LOG_OPEN);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo;

     FILE *fp_emerged_halos;
     char  filename_emerged_halos_out[256];
     sprintf(filename_out_dir,"%s/emerged_halos",filename_out_root_dir);
     if(i_snap==0)
        mkdir(filename_out_dir,02755);
     sprintf(filename_emerged_halos_out,"%s/%03d.txt",filename_out_dir,i_snap);
     fp_emerged_halos=fopen(filename_emerged_halos_out,"w");

     FILE *fp_simple_halos;
     char  filename_simple_halos_out[256];
     sprintf(filename_out_dir,"%s/simple_halos",filename_out_root_dir);
     if(i_snap==0)
        mkdir(filename_out_dir,02755);
     sprintf(filename_simple_halos_out,"%s/%03d.txt",filename_out_dir,i_snap);
     fp_simple_halos=fopen(filename_simple_halos_out,"w");

     // Loop over each halo for this snapshot ...
     current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){
        tree_node_info *main_progenitor;
        tree_node_info *merging_progenitor;

        // Set the main progenitor ...
        main_progenitor=current_halo->progenitor_first;
        if(main_progenitor!=NULL){
           int i_np_0=0; 
           if(main_progenitor->n_particles<bin_np_lo[i_np_0])
              i_np_0=-1;
           else{
              while(main_progenitor->n_particles>=bin_np_hi[i_np_0] && i_np_0<n_np) 
                 i_np_0++;
           }

           // Loop over each progenitor that is merging with this main progenitor ...
           if(i_np_0>=0 && i_np_0<n_np){
              // Check if the main progenitor is a fragmented halo
              int flag_fragmented=check_mode_for_flag(main_progenitor->tree_case,TREE_CASE_FRAGMENTED_RETURNED) ||
                                  check_mode_for_flag(main_progenitor->tree_case,TREE_CASE_FRAGMENTED_STRAYED)||
                                  check_mode_for_flag(main_progenitor->tree_case,TREE_CASE_FRAGMENTED_EXCHANGED);

              // Populate counters
              int i_bin=i_np_0*n_np+i_np_0;
              if(flag_fragmented)
                 n_fragmented[i_snap][i_bin]++;

              merging_progenitor=main_progenitor->progenitor_next;
              while(merging_progenitor!=NULL){
                 int i_np_i=0;
                 if(merging_progenitor->n_particles<bin_np_lo[i_np_i])
                    i_np_i=-1;
                 else{
                    while(merging_progenitor->n_particles>=bin_np_hi[i_np_i] && i_np_i<n_np) 
                       i_np_i++;
                 }
                 // Process each merger ...
                 if(i_np_i>=0 && i_np_i<n_np){
                    // Find the start of this merging halo's main progenitor line
                    tree_node_info *current;
                    int i_snap_start;
                    int di;
                    current=merging_progenitor;
                    while(current!=NULL){
                       i_snap_start=current->snap_tree;
                       current=current->progenitor_first;
                    }
                    di=i_snap-i_snap_start-1;

                    // Determine if this halo is actually a fragmented halo
                    int flag_fragmented=check_mode_for_flag(merging_progenitor->tree_case,TREE_CASE_FRAGMENTED_RETURNED) ||
                                        check_mode_for_flag(merging_progenitor->tree_case,TREE_CASE_FRAGMENTED_STRAYED)||
                                        check_mode_for_flag(merging_progenitor->tree_case,TREE_CASE_FRAGMENTED_EXCHANGED);

                    // Populate counters
                    int i_bin=i_np_0*n_np+i_np_i;
                    if(flag_fragmented){
                       n_fragmented[i_snap][i_bin]++;
                       n_fragmented_ages[i_snap][i_bin*i_snap+di]++;
                    }
                    else{
                       n_mergers[i_snap][i_bin]++;
                       n_mergers_ages[i_snap][i_bin*i_snap+di]++;
                    }
                 }
                 merging_progenitor=merging_progenitor->progenitor_next;
              }
           }
        }

        // Perform simpler counts related to this snapshot's halos
        int i_np=0;
        if(current_halo->n_particles<bin_np_lo[i_np])
           i_np=-1;
        else{
           while(current_halo->n_particles>=bin_np_hi[i_np] && i_np<n_np) 
              i_np++;
        }
        if(i_np>=0 && i_np<n_np){
           // Count dropped halos
           int flag_dropped=check_mode_for_flag(current_halo->tree_case,TREE_CASE_DROPPED);
           if(flag_dropped){
              int di=current_halo->descendant->snap_tree-current_halo->snap_tree-1;
              n_dropped[i_snap][i_np]++;
              n_dropped_ages[i_snap][i_np*i_snap+di]++;
           }
           // Count emerged halos
           int flag_emerged=FALSE;
           if(current_halo->progenitor_first!=NULL){
              if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_EMERGED))
                 flag_emerged=TRUE;
           }
           if(flag_emerged){
              int di=current_halo->snap_tree-current_halo->progenitor_first->snap_tree-1;
              n_emerged[i_snap][i_np]++;
              n_emerged_ages[i_snap][i_np*i_snap+di]++;
              fprintf(fp_emerged_halos,"%d %d %d %d %d %d %e\n",
                                       trees->snap_list[current_halo->progenitor_first->snap_tree],
                                       current_halo->progenitor_first->file_index,
                                       trees->snap_list[i_snap],
                                       current_halo->file_index,
                                       current_halo->progenitor_first->n_particles,
                                       current_halo->n_particles,
                                       match_score_local[current_halo->progenitor_first->snap_tree][current_halo->progenitor_first->neighbour_index]);
           }
           else if(check_mode_for_flag(current_halo->tree_case,TREE_CASE_SIMPLE) && current_halo->descendant!=NULL){
              fprintf(fp_simple_halos,"%d %d %d %d %d %d %e\n",
                                      trees->snap_list[i_snap],
                                      current_halo->file_index,
                                      trees->snap_list[current_halo->descendant->snap_tree],
                                      current_halo->descendant->file_index,
                                      current_halo->n_particles,
                                      current_halo->descendant->n_particles,
                                      match_score_local[current_halo->snap_tree][current_halo->neighbour_index]);
           }
           // Count the total number of halos
           n_halos[i_snap][i_np]++;
        }

        current_halo=current_halo->next_neighbour;
     }
     SID_Allreduce(SID_IN_PLACE,n_fragmented[i_snap],n_bins,SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,n_mergers[i_snap],   n_bins,SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,n_halos[i_snap],     n_np,  SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,n_emerged[i_snap],   n_np,  SID_INT,SID_SUM,SID.COMM_WORLD);
     SID_Allreduce(SID_IN_PLACE,n_dropped[i_snap],   n_np,  SID_INT,SID_SUM,SID.COMM_WORLD);
     if(i_snap>0){
        SID_Allreduce(SID_IN_PLACE,n_fragmented_ages[i_snap],i_snap*n_bins,SID_INT,SID_SUM,SID.COMM_WORLD);
        SID_Allreduce(SID_IN_PLACE,n_mergers_ages[i_snap],   i_snap*n_bins,SID_INT,SID_SUM,SID.COMM_WORLD);
        SID_Allreduce(SID_IN_PLACE,n_emerged_ages[i_snap],   i_snap*n_np,  SID_INT,SID_SUM,SID.COMM_WORLD);
        SID_Allreduce(SID_IN_PLACE,n_dropped_ages[i_snap],   i_snap*n_np,  SID_INT,SID_SUM,SID.COMM_WORLD);
     }
     fclose(fp_emerged_halos);
     fclose(fp_simple_halos);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Write merger counts
  sprintf(filename_out_dir,"%s/mergers",filename_out_root_dir);
  mkdir(filename_out_dir,02755);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
    FILE *fp;
    char  filename[256];
    sprintf(filename,"%s/%03d.txt",filename_out_dir,i_snap);
    fp=fopen(filename,"w");
    int i_bin;
    for(int i_np_0=0,i_bin=0;i_np_0<n_np;i_np_0++){
       for(int i_np_i=0;i_np_i<n_np;i_np_i++,i_bin++){
          if(i_np_0<i_np_i)
             fprintf(fp,"%d %d %d\n",
                        bin_np_lo[i_np_0],
                        bin_np_lo[i_np_i],
                        n_mergers[i_snap][i_np_i*n_np+i_np_0]);
          else
             fprintf(fp,"%d %d %d\n",
                        bin_np_lo[i_np_0],
                        bin_np_lo[i_np_i],
                        n_mergers[i_snap][i_np_0*n_np+i_np_i]);
       }
       fprintf(fp,"\n");
    }
    fclose(fp);
  }

  // Write fragmented halo counts
  sprintf(filename_out_dir,"%s/fragmented",filename_out_root_dir);
  mkdir(filename_out_dir,02755);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
    FILE *fp;
    char  filename[256];
    sprintf(filename,"%s/%03d.txt",filename_out_dir,i_snap);
    fp=fopen(filename,"w");
    int i_bin;
    for(int i_np_0=0,i_bin=0;i_np_0<n_np;i_np_0++){
       for(int i_np_i=0;i_np_i<n_np;i_np_i++,i_bin++){
          if(i_np_0<i_np_i)
             fprintf(fp,"%d %d %d\n",
                        bin_np_lo[i_np_0],
                        bin_np_lo[i_np_i],
                        n_fragmented[i_snap][i_np_i*n_np+i_np_0]);
          else
             fprintf(fp,"%d %d %d\n",
                        bin_np_lo[i_np_0],
                        bin_np_lo[i_np_i],
                        n_fragmented[i_snap][i_np_0*n_np+i_np_i]);
       }
       fprintf(fp,"\n");
    }
    fclose(fp);
  }

  // Write collapsed halo counts
  SID_log("Writing results...",SID_LOG_OPEN);
  int collect_size;
  for(collect_size=1;collect_size<=3;collect_size++){
     char  filename[256];
     FILE *fp_z;
     FILE *fp_t;
     int   collect_start=0;
     int   n_collect=0;
     int   i_collect;
     int   i_column;
     int   j_np;
     for(j_np=collect_start;j_np<n_np;j_np+=collect_size) n_collect++;
     int *n_dropped_z   =(int *)SID_malloc(sizeof(int)*n_collect);
     int *n_emerged_z   =(int *)SID_malloc(sizeof(int)*n_collect);
     int *n_fragmented_z=(int *)SID_malloc(sizeof(int)*n_collect);
     int *n_mergers_z   =(int *)SID_malloc(sizeof(int)*n_collect);
     int *n_halos_z     =(int *)SID_malloc(sizeof(int)*n_collect);
     sprintf(filename_out_dir,"%s/ages",filename_out_root_dir);
     mkdir(filename_out_dir,02755);
     for(i_snap=1;i_snap<trees->n_snaps;i_snap++){
        sprintf(filename,"%s/%03d_%03d.txt",filename_out_dir,collect_size,i_snap);
        fp_t=fopen(filename,"w");
        fprintf(fp_t,"# Halo age distributions for {%s} at z=%le\n",filename_out_root,trees->z_list[i_snap]);
        fprintf(fp_t,"#\n");
        i_column=1;
        fprintf(fp_t,"# Column (%02d): Delta [Gyrs]\n",i_column++);
        for(j_np=collect_start;j_np<n_np;j_np+=collect_size){
           int i_lo=j_np;
           int i_hi=MIN(j_np+collect_size-1,n_np-1);
           fprintf(fp_t,"#        (%02d): n_dropped    [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
           fprintf(fp_t,"#        (%02d): n_emerged    [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
           fprintf(fp_t,"#        (%02d): n_fragmented [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
           fprintf(fp_t,"#        (%02d): n_mergers    [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
        }
        int j_snap;
        for(j_snap=0;j_snap<i_snap;j_snap++){
          fprintf(fp_t,"%le",(trees->t_list[i_snap]-trees->t_list[i_snap-j_snap-1])/(1e9*S_PER_YEAR));
          for(j_np=collect_start;j_np<n_np;){
             int n_dropped_collect;
             int n_emerged_collect;
             int n_fragmented_collect;
             int n_mergers_collect;
             n_dropped_collect   =0;
             n_emerged_collect   =0;
             n_fragmented_collect=0;
             n_mergers_collect   =0;
             for(int j_collect=0;j_collect<collect_size && j_np<n_np;j_collect++,j_np++){
                n_emerged_collect+=n_emerged_ages[i_snap][j_np*i_snap+j_snap];
                n_dropped_collect+=n_dropped_ages[i_snap][j_np*i_snap+j_snap];
                for(i_np=0;i_np<n_np;i_np++){
                   n_fragmented_collect+=n_fragmented_ages[i_snap][(i_np*n_np+j_np)*i_snap+j_snap];
                   n_mergers_collect   +=n_mergers_ages[i_snap][(i_np*n_np+j_np)*i_snap+j_snap];
                }
             }
             fprintf(fp_t," %5d %5d %5d %5d",n_dropped_collect,n_emerged_collect,n_fragmented_collect,n_mergers_collect);
          }
          fprintf(fp_t,"\n");
        }
        fclose(fp_t);
     }

     sprintf(filename,"%s/%03d_z.txt",filename_out_root_dir,collect_size);
     fp_z=fopen(filename,"w");
     fprintf(fp_z,"# Halo counts for {%s}\n",filename_out_root);
     fprintf(fp_z,"#\n");
     i_column=1;
     fprintf(fp_z,"# Column (%02d): redshift\n",i_column++);
     fprintf(fp_z,"#        (%02d): t [Gyrs]\n",i_column++);
     for(int j_np=collect_start;j_np<n_np;j_np+=collect_size){
        int i_lo=j_np;
        int i_hi=MIN(j_np+collect_size-1,n_np-1);
        fprintf(fp_z,"#        (%02d): n_dropped    [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
        fprintf(fp_z,"#        (%02d): n_emerged    [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
        fprintf(fp_z,"#        (%02d): n_fragmented [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
        fprintf(fp_z,"#        (%02d): n_mergers    [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
        fprintf(fp_z,"#        (%02d): n_halos      [%d<=n_p_i<%d]\n",i_column++,bin_np_lo[i_lo],bin_np_hi[i_hi]);
     }
     sprintf(filename_out_dir,"%s/ni",filename_out_root_dir);
     mkdir(filename_out_dir,02755);
     for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
       fprintf(fp_z,"%le %le",trees->z_list[i_snap],trees->t_list[i_snap]/(1e9*S_PER_YEAR));
       FILE *fp_ni;
       sprintf(filename,"%s/%03d_%03d.txt",filename_out_dir,collect_size,i_snap);
       fp_ni=fopen(filename,"w");
       fprintf(fp_ni,"# Halo counts for {%s}, z=%le\n",filename_out_root,trees->z_list[i_snap]);
       fprintf(fp_ni,"#\n");
       int i_column=1;
       fprintf(fp_ni,"# Column (%02d): n_p bin (lo)\n",       i_column++);
       fprintf(fp_ni,"#        (%02d): n_p bin (hi)\n",       i_column++);
       fprintf(fp_ni,"#        (%02d): n_halos\n",            i_column++);
       fprintf(fp_ni,"#        (%02d): n_dropped\n",          i_column++);
       fprintf(fp_ni,"#        (%02d): n_emerged\n",          i_column++);
       fprintf(fp_ni,"#        (%02d): n_fragmented(n_p_i)\n",i_column++);
       fprintf(fp_ni,"#        (%02d): n_mergers(n_p_i)\n",   i_column++);
       int j_np;
       int n_fragmented_collect=0;
       int n_mergers_collect   =0;
       for(int i_np=0;i_np<n_np;i_np++){
          fprintf(fp_ni,"%7d %7d %6d %6d %6d",bin_np_lo[i_np],bin_np_hi[i_np],n_halos[i_snap][i_np],n_dropped[i_snap][i_np],n_emerged[i_snap][i_np]);
          n_fragmented_collect=0;
          n_mergers_collect   =0;
          for(j_np=0;j_np<n_np;j_np++){
             n_fragmented_collect+=n_fragmented[i_snap][j_np*n_np+i_np];
             n_mergers_collect   +=n_mergers[i_snap][j_np*n_np+i_np];
          }
          fprintf(fp_ni," %5d %5d\n",n_fragmented_collect,n_mergers_collect);
       }
       fclose(fp_ni);
       for(i_collect=0;i_collect<n_collect;i_collect++){
          n_dropped_z[i_collect]   =0;
          n_emerged_z[i_collect]   =0;
          n_fragmented_z[i_collect]=0;
          n_mergers_z[i_collect]   =0;
          n_halos_z[i_collect]     =0;
       }
       for(int i_np=0;i_np<n_np;i_np++){
          for(j_np=collect_start,i_collect=0;j_np<n_np;i_collect++){
             for(int j_collect=0;j_collect<collect_size && j_np<n_np;j_collect++,j_np++){
                n_fragmented_z[i_collect]+=n_fragmented[i_snap][i_np*n_np+j_np];
                n_mergers_z[i_collect]   +=n_mergers[i_snap][i_np*n_np+j_np];
             }
          }
       }
       for(j_np=collect_start,i_collect=0;j_np<n_np;i_collect++){
          for(int j_collect=0;j_collect<collect_size && j_np<n_np;j_collect++,j_np++){
             n_dropped_z[i_collect]+=n_dropped[i_snap][j_np];
             n_emerged_z[i_collect]+=n_emerged[i_snap][j_np];
             n_halos_z[i_collect]  +=n_halos[i_snap][j_np];
          }
       }
       for(i_collect=0;i_collect<n_collect;i_collect++)
          fprintf(fp_z," %5d %5d %5d %5d %5d",n_dropped_z[i_collect],n_emerged_z[i_collect],n_fragmented_z[i_collect],n_mergers_z[i_collect],n_halos_z[i_collect]);
       fprintf(fp_z,"\n");
     }
     fclose(fp_z);
     SID_free(SID_FARG n_dropped_z);
     SID_free(SID_FARG n_emerged_z);
     SID_free(SID_FARG n_fragmented_z);
     SID_free(SID_FARG n_mergers_z);
     SID_free(SID_FARG n_halos_z);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  SID_free(SID_FARG bin_np_lo);
  SID_free(SID_FARG bin_np_hi);
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     SID_free(SID_FARG n_halos[i_snap]);
     SID_free(SID_FARG n_dropped[i_snap]);
     SID_free(SID_FARG n_emerged[i_snap]);
     SID_free(SID_FARG n_fragmented[i_snap]);
     SID_free(SID_FARG n_mergers[i_snap]);
     SID_free(SID_FARG n_dropped_ages[i_snap]);
     SID_free(SID_FARG n_emerged_ages[i_snap]);
     SID_free(SID_FARG n_fragmented_ages[i_snap]);
     SID_free(SID_FARG n_mergers_ages[i_snap]);
  }
  SID_free(SID_FARG n_halos);
  SID_free(SID_FARG n_dropped);
  SID_free(SID_FARG n_emerged);
  SID_free(SID_FARG n_fragmented);
  SID_free(SID_FARG n_mergers);
  SID_free(SID_FARG n_dropped_ages);
  SID_free(SID_FARG n_emerged_ages);
  SID_free(SID_FARG n_fragmented_ages);
  SID_free(SID_FARG n_mergers_ages);

  SID_log("Done.",SID_LOG_CLOSE);
}

