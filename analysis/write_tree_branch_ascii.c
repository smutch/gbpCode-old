#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void write_tree_branch_ascii(tree_info *trees,tree_node_info *halo,const char *filename_out,const char *trees_name){

  SID_log("Writing branch to file {%s}...",SID_LOG_OPEN,filename_out);
  if(halo!=NULL){
     // Find branch markers
     tree_markers_info markers;
     find_treenode_markers(trees,halo,&markers);

     // Perform write
     int   i_column=1;
     FILE *fp_out  =fopen(filename_out,"w");
     fprintf(fp_out,"# Main progenitor line for halo ID No. %d of {%s}\n",halo->halo_ID,trees_name);
     fprintf(fp_out,"#\n");
     fprintf(fp_out,"# Important halo markers (snapshot,snap_tree,halo_index):\n");
     fprintf(fp_out,"#    selected halo        = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,halo),
                    fetch_treenode_snap_tree(trees,halo),
                    fetch_treenode_file_index(trees,halo));
     fprintf(fp_out,"#    main progenitor      = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.main_progenitor),
                    fetch_treenode_snap_tree(trees,markers.main_progenitor),
                    fetch_treenode_file_index(trees,markers.main_progenitor));
     fprintf(fp_out,"#    branch root          = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.branch_root),
                    fetch_treenode_snap_tree(trees,markers.branch_root),
                    fetch_treenode_file_index(trees,markers.branch_root));
     fprintf(fp_out,"#    branch leaf          = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.branch_leaf),
                    fetch_treenode_snap_tree(trees,markers.branch_leaf),
                    fetch_treenode_file_index(trees,markers.branch_leaf));
     fprintf(fp_out,"#    descendant           = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.descendant),
                    fetch_treenode_snap_tree(trees,markers.descendant),
                    fetch_treenode_file_index(trees,markers.descendant));
     fprintf(fp_out,"#    1st became satellite = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.first_became_satellite),
                    fetch_treenode_snap_tree(trees,markers.first_became_satellite),
                    fetch_treenode_file_index(trees,markers.first_became_satellite));
     fprintf(fp_out,"#    joined current group = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.joined_current_parent),
                    fetch_treenode_snap_tree(trees,markers.joined_current_parent),
                    fetch_treenode_file_index(trees,markers.joined_current_parent));
     fprintf(fp_out,"#    peak mass            = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.peak_mass),
                    fetch_treenode_snap_tree(trees,markers.peak_mass),
                    fetch_treenode_file_index(trees,markers.peak_mass));
     fprintf(fp_out,"#    half-peak-mass       = (%d,%d,%d)\n",
                    fetch_treenode_snapshot(trees,markers.half_peak_mass),
                    fetch_treenode_snap_tree(trees,markers.half_peak_mass),
                    fetch_treenode_file_index(trees,markers.half_peak_mass));
     fprintf(fp_out,"#\n");
     fprintf(fp_out,"# Column (%02d): snapshot\n",          i_column++);
     fprintf(fp_out,"#        (%02d): snapshot index\n",    i_column++);
     fprintf(fp_out,"#        (%02d): halo index\n",        i_column++);
     fprintf(fp_out,"#        (%02d): halo ID\n",           i_column++);
     fprintf(fp_out,"#        (%02d): group ID\n",          i_column++);
     fprintf(fp_out,"#        (%02d): central ID\n",        i_column++);
     fprintf(fp_out,"#        (%02d): M_vir [h^-1 M_sol]\n",i_column++);
     fprintf(fp_out,"#\n");
     tree_node_info *current_halo=markers.branch_root;
     while(current_halo!=NULL){
        double M_vir=fetch_treenode_Mvir(trees,current_halo);
        fprintf(fp_out,"%3d %3d %7d %7d %7d %7d %le\n",
                       fetch_treenode_snapshot(trees,current_halo),
                       fetch_treenode_snap_tree(trees,current_halo),
                       fetch_treenode_file_index(trees,current_halo),
                       current_halo->halo_ID,
                       current_halo->parent->halo_ID,
                       current_halo->parent->substructure_first->halo_ID,
                       M_vir);
        current_halo=current_halo->progenitor_first;
     }
     fclose(fp_out);
  }
  else
     SID_log("skipping NULL halo...",SID_LOG_CONTINUE);
  SID_log("Done.",SID_LOG_CLOSE);
}

