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
#include <assert.h>

void compute_trees_emerged_halo_analysis(tree_info *trees,char *filename_out_root){

  // Compute merger rates ...
  SID_log("Performing emerged halo analysis...",SID_LOG_OPEN|SID_LOG_TIMER);
  int   i_snap;
  int   i_np;
  int   i_xi;

  // Set output directory
  FILE *fp_emerged_sums;
  FILE *fp_emerged_list;
  char  filename_emerged_sums_out[256];
  char  filename_emerged_list_out[256];
  char  filename_out_root_dir[MAX_FILENAME_LENGTH];
  char  filename_out_dir[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root_dir,    "%s_tree_analysis/",filename_out_root);
  sprintf(filename_out_dir,         "%s/emerged_halos", filename_out_root_dir);
  sprintf(filename_emerged_sums_out,"%s/sums.txt", filename_out_dir);
  sprintf(filename_emerged_list_out,"%s/list.txt", filename_out_dir);
  mkdir(filename_out_root_dir,02755);
  mkdir(filename_out_dir,     02755);
  fp_emerged_sums=fopen(filename_emerged_sums_out,"w");
  fp_emerged_list=fopen(filename_emerged_list_out,"w");

  // Allocate a couple temporary arrays
  tree_markers_info  *markers_list=(tree_markers_info  *)SID_malloc(sizeof(tree_markers_info)*trees->n_snaps);
  tree_node_info    **emerged_list=(tree_node_info    **)SID_malloc(sizeof(tree_node_info *) *trees->n_snaps);

  // Loop over each snapshot ...
  for(i_snap=0;i_snap<trees->n_snaps;i_snap++){
     tree_node_info *current_halo=trees->first_neighbour_subgroups[i_snap];
     while(current_halo!=NULL){

        // Process each new branch
        if(check_treenode_if_branch_start(current_halo)){
           // Scan the branch
           int n_emerge=0;
           tree_node_info *current_progenitor=current_halo;
           while(current_progenitor!=NULL){
              // Process emerged halos in this branch's main progenitor line
              if(check_mode_for_flag(current_progenitor->tree_case,TREE_CASE_EMERGED)){
                 emerged_list[n_emerge]=current_halo;
                 find_treenode_markers(trees,current_progenitor,&(markers_list[n_emerge]));
                 n_emerge++;
              }
              current_progenitor=current_progenitor->progenitor_first;
           }
           // this branch has emerged halos in it.  Print statistics.
           if(n_emerge>0){
              fprintf(fp_emerged_sums,"%le %le %le %d %d %d\n",
                                      trees->z_list[current_halo->snap_tree],
                                      fetch_treenode_Mvir(trees,current_halo),
                                      fetch_treenode_Mvir(trees,current_halo->parent),
                                      fetch_treenode_n_particles (trees,current_halo),
                                      fetch_treenode_n_particles (trees,current_halo->parent),
                                      n_emerge);
              for(int i_emerge=0;i_emerge<n_emerge;i_emerge++)
                 fprintf(fp_emerged_list,"%d %le %le %le %le %le %le %le %d %d %d %d %d %d %d %le %le %le %le %le %le %le %le %le\n",
                                         i_emerge,
                                         fetch_treenode_Mvir   (trees,emerged_list[i_emerge]),
                                         fetch_treenode_Mvir   (trees,emerged_list[i_emerge]->parent),
                                         fetch_treenode_Mvir   (trees,markers_list[i_emerge].progenitor_main),
                                         fetch_treenode_Mvir   (trees,markers_list[i_emerge].descendant_last),
                                         fetch_treenode_Mvir   (trees,markers_list[i_emerge].progenitor_accretion_first),
                                         fetch_treenode_Mvir   (trees,markers_list[i_emerge].progenitor_accretion_last),
                                         fetch_treenode_Mvir   (trees,markers_list[i_emerge].progenitor_formation),
                                         fetch_treenode_n_particles    (trees,emerged_list[i_emerge]),
                                         fetch_treenode_n_particles    (trees,emerged_list[i_emerge]->parent),
                                         fetch_treenode_n_particles    (trees,markers_list[i_emerge].progenitor_main),
                                         fetch_treenode_n_particles    (trees,markers_list[i_emerge].descendant_last),
                                         fetch_treenode_n_particles    (trees,markers_list[i_emerge].progenitor_accretion_first),
                                         fetch_treenode_n_particles    (trees,markers_list[i_emerge].progenitor_accretion_last),
                                         fetch_treenode_n_particles    (trees,markers_list[i_emerge].progenitor_formation),
                                         fetch_treenode_z      (trees,emerged_list[i_emerge]),
                                         fetch_treenode_z      (trees,markers_list[i_emerge].descendant_last),
                                         fetch_treenode_z      (trees,markers_list[i_emerge].progenitor_accretion_first),
                                         fetch_treenode_z      (trees,markers_list[i_emerge].progenitor_accretion_last),
                                         fetch_treenode_z      (trees,markers_list[i_emerge].progenitor_formation),
                                         fetch_treenode_delta_t(trees,current_halo,markers_list[i_emerge].descendant_last)/S_PER_YEAR,
                                         fetch_treenode_delta_t(trees,current_halo,markers_list[i_emerge].progenitor_accretion_first)/S_PER_YEAR,
                                         fetch_treenode_delta_t(trees,current_halo,markers_list[i_emerge].progenitor_accretion_last)/S_PER_YEAR,
                                         fetch_treenode_delta_t(trees,current_halo,markers_list[i_emerge].progenitor_formation)/S_PER_YEAR);
           }
        }
        current_halo=current_halo->next_neighbour;
     }
  }
  fclose(fp_emerged_sums);
  fclose(fp_emerged_list);

  // Clean-up
  SID_free(SID_FARG emerged_list);
  SID_free(SID_FARG markers_list);

  SID_log("Done.",SID_LOG_CLOSE);
}

