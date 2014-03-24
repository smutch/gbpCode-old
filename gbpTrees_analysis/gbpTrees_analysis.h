#ifndef ZFIRE_AWAKE
#define ZFIRE_AWAKE

// V Preprocessor definitions V
// A Preprocessor definitions A

// V --- Datatype definitions --- V
typedef struct tree_markers_info tree_markers_info;
struct tree_markers_info{
  int             flag_halo_is_main_progenitor;
  tree_node_info *leaf;
  tree_node_info *last_instance;
  tree_node_info *progenitor_main;
  tree_node_info *progenitor_accretion_first;
  tree_node_info *progenitor_accretion_last;
  tree_node_info *progenitor_formation;
};

// A --- Datatype definitions --- A

#ifdef __cplusplus
extern "C" {
#endif   
// V --- ANSI-C function definitions --- V
void ZFIRE_compute_catalog(tree_info *trees,const char *catalog_name);
void ZFIRE_compute_catalog_main(tree_info *trees,const char *catalog_name);
int check_treenode_if_main_progenitor(tree_node_info *halo);
halo_properties_info *fetch_treenode_properties(tree_info *trees,tree_node_info *halo);
int find_treesnap_z(tree_info *trees,double z_obs_exact);
int find_treenode_main_progenitor(tree_info *trees,tree_node_info *halo,tree_node_info **main_progenitor);
int find_treenode_leaf(tree_info *trees,tree_node_info *halo,tree_node_info **leaf);
int find_treenode_last_instance(tree_info *trees,tree_node_info *halo,tree_node_info **last_descendant);
int find_treenode_formation(tree_info *trees,tree_node_info *halo,double f,tree_node_info **formation_progenitor);
int find_treenode_accretion(tree_info       *trees,
                            tree_node_info  *halo,
                            tree_node_info **first_accretion_progenitor,
                            tree_node_info **last_accretion_progenitor);
int find_treenode_markers(tree_info *trees,tree_node_info *halo,tree_markers_info *markers);
int fetch_treenode_snapshot(tree_info *trees,tree_node_info *halo);
int fetch_treenode_snap_tree(tree_info *trees,tree_node_info *halo);
int fetch_treenode_file_index(tree_info *trees,tree_node_info *halo);
double fetch_treenode_Mvir(tree_info *trees,tree_node_info *halo);
void write_tree_tracks(tree_info *trees,tree_node_info **list_in,int n_list_in,int mode,const char *catalog_name);
void average_tree_tracks(const char *catalog_name);

// A --- ANSI-C function definitions --- A
#ifdef __cplusplus
}
#endif

#endif
