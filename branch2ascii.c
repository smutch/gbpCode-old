#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>
#include <gbpZFIRE.h>

int main(int argc, char *argv[]){

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  char filename_SSimPL_dir[MAX_FILENAME_LENGTH];
  char filename_halo_version_root[MAX_FILENAME_LENGTH];
  char filename_trees_root[MAX_FILENAME_LENGTH];
  char filename_trees_name[MAX_FILENAME_LENGTH];
  char filename_halos_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_dir,       argv[1]);
  strcpy(filename_halo_version_root,argv[2]);
  strcpy(filename_trees_name,       argv[3]);
  int halo_ID                 =atoi(argv[4]);

  sprintf(filename_trees_root,"%s/trees/%s",filename_SSimPL_dir,filename_trees_name);
  sprintf(filename_halos_root,"%s/halos/%s",filename_SSimPL_dir,filename_halo_version_root);

  SID_log("Writing an ascii version of halo_id=%d's main progenitor branch...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Perform analysis
  tree_info *trees;
  read_trees(filename_trees_root,
             filename_halos_root,
             TREE_MODE_DEFAULT,
             &trees);

  // Read ancillary data
  read_trees_catalogs(trees,
                      filename_SSimPL_dir,
                      filename_halo_version_root,
                      READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS);

  // Select forest corresponding to given halo id
  tree_node_info *first_instance=NULL;
  //for(int i_forest=0;i_forest<trees->n_forests_local && first_instance==NULL;i_forest++){
  for(int i_forest=(trees->n_forests_local-1);i_forest>=0 && first_instance==NULL;i_forest--){
     tree_node_info *current_halo=trees->first_in_forest_subgroups[i_forest];
     while(current_halo && first_instance==NULL){
        if(current_halo->halo_ID==halo_ID)
           first_instance=current_halo;
        current_halo=current_halo->next_in_forest;
     }
  }

  // Create filename
  char  filename_out[MAX_FILENAME_LENGTH];
  sprintf(filename_out,"halo_id_%09d.ascii",halo_ID);
  SID_log("Writing results to {%s}...",SID_LOG_OPEN,filename_out);

  // Find branch markers
  tree_markers_info markers;
  find_treenode_markers(trees,first_instance,       &markers);
  find_treenode_markers(trees,markers.last_instance,&markers);

  // Perform write
  FILE *fp_out;
  int   i_column=1;
  fp_out=fopen(filename_out,"w");
  fprintf(fp_out,"# Main progenitor line for halo ID No. %d of {%s}\n",halo_ID,filename_trees_name);
  fprintf(fp_out,"#\n");
/*
  fprintf(fp_out,"# Important halos (snapshot,snap_tree,halo_index):\n");
  fprintf(fp_out,"#    last instance              = (%d,%d,%d)\n",
                 fetch_treenode_snapshot(trees,markers.last_instance),
                 fetch_treenode_snap_tree(trees,markers.last_instance),
                 fetch_treenode_file_index(trees,markers.last_instance));
  fprintf(fp_out,"#    descendant                 = (%d,%d,%d)\n",
                 fetch_treenode_snapshot(trees,markers.last_instance->descendant),
                 fetch_treenode_snap_tree(trees,markers.last_instance->descendant),
                 fetch_treenode_file_index(trees,markers.last_instance->descendant));
  fprintf(fp_out,"#    progenitor_accretion_first = (%d,%d,%d)\n",
                 fetch_treenode_snapshot(trees,markers.progenitor_accretion_first),
                 fetch_treenode_snap_tree(trees,markers.progenitor_accretion_first),
                 fetch_treenode_file_index(trees,markers.progenitor_accretion_first));
  fprintf(fp_out,"#    progenitor_accretion_last  = (%d,%d,%d)\n",
                 fetch_treenode_snapshot(trees,markers.progenitor_accretion_last),
                 fetch_treenode_snap_tree(trees,markers.progenitor_accretion_last),
                 fetch_treenode_file_index(trees,markers.progenitor_accretion_last));
  fprintf(fp_out,"#    progenitor_formation       = (%d,%d,%d)\n",
                 fetch_treenode_snapshot(trees,markers.progenitor_formation),
                 fetch_treenode_snap_tree(trees,markers.progenitor_formation),
                 fetch_treenode_file_index(trees,markers.progenitor_formation));
  fprintf(fp_out,"#\n");
*/
  fprintf(fp_out,"# Column (%02d): snapshot\n",          i_column++);
  fprintf(fp_out,"#        (%02d): snap_tree\n",         i_column++);
  fprintf(fp_out,"#        (%02d): file index\n",        i_column++);
  fprintf(fp_out,"#        (%02d): halo ID\n",           i_column++);
  fprintf(fp_out,"#        (%02d): group ID\n",          i_column++);
  fprintf(fp_out,"#        (%02d): central ID\n",        i_column++);
  fprintf(fp_out,"#        (%02d): M_vir [h^-1 M_sol]\n",i_column++);
  fprintf(fp_out,"#\n");
  tree_node_info *current_halo=markers.last_instance;
  while(current_halo!=NULL){
     fprintf(fp_out,"%3d %3d %7d %7d %7d %7d %le\n",
                    fetch_treenode_snapshot(trees,current_halo),
                    fetch_treenode_snap_tree(trees,current_halo),
                    fetch_treenode_file_index(trees,current_halo),
                    current_halo->halo_ID,
                    current_halo->parent->halo_ID,
                    current_halo->parent->substructure_first->halo_ID,
                    fetch_treenode_Mvir(trees,current_halo));
     current_halo=current_halo->progenitor_first;
  }
  fclose(fp_out);
  SID_log("Done.",SID_LOG_CLOSE);

  // Clean-up
  free_trees(&trees);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

