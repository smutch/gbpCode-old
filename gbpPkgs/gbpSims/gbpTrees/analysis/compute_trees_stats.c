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

void compute_trees_stats(tree_info *trees){

  // Make sure the output directory exists
  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s/treenode",trees->filename_root_analysis);
  mkdir(trees->filename_root_analysis,02755);

  // Loop twice; once for groups and once for subgroups
  for(int i_type=0;i_type<2;i_type++){
     //// Precompute marker statistics
     //if(i_type==0)
     //   precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     //else
     //   precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);

     // Fetch the array that lists the first neighbour of each snapshot
     tree_node_info **neighbour_list_start=NULL;
     char group_prefix[8];
     if(i_type==0){
        sprintf(group_prefix,"");
        neighbour_list_start=trees->first_neighbour_groups;
     }
     else{
        sprintf(group_prefix,"sub");
        neighbour_list_start=trees->first_neighbour_subgroups;
     }
   
     // Count the number of mergers
     SID_log("Counting mergers...",SID_LOG_OPEN|SID_LOG_TIMER);
     int sum_mergers =0;
     int sum_prog_nth=0;
     for(int i_snap=0;i_snap<trees->n_snaps;i_snap++){
        tree_node_info *current_halo=neighbour_list_start[i_snap];
        int n_mergers =0;
        int n_prog_nth=0;
        while(current_halo!=NULL){
           // Process each new merger
           if(check_treenode_if_merger(current_halo))
              n_mergers++;
           // Count nth progenitors
           int i_prog=0;
           tree_node_info *current_prog=current_halo->progenitor_first;
           while(current_prog!=NULL){
              if(i_prog!=0)
                 n_prog_nth++;
              i_prog++;
              current_prog=current_prog->progenitor_next;
           }
           current_halo=current_halo->next_neighbour;
        }
        sum_mergers +=n_mergers;
        sum_prog_nth+=n_prog_nth;
        printf("%3d %3d %3d\n",trees->snap_list[i_snap],n_mergers,n_prog_nth);
     }
     printf("Sum mergers  =%d\n",sum_mergers);
     printf("Sum nth progs=%d\n",sum_prog_nth);
     SID_log("Done.",SID_LOG_CLOSE);     
   
     //// Free catalogs and precomputed markers
     //if(i_type==0)
     //   free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     //else
     //   free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
  }
}

