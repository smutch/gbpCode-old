#ifndef GBPTREES_ANALYSIS_AWAKE
#define GBPTREES_ANALYSIS_AWAKE

// V Preprocessor definitions V
// A Preprocessor definitions A

// V --- Datatype definitions --- V
typedef struct tree_markers_info tree_markers_info;
struct tree_markers_info{
  int             flag_halo_is_main_progenitor;
  tree_node_info *branch_leaf;
  tree_node_info *branch_root;
  tree_node_info *descendant;
  tree_node_info *main_progenitor;
  tree_node_info *first_became_satellite;
  tree_node_info *joined_current_parent;
  tree_node_info *peak_mass;
  tree_node_info *half_peak_mass;
};
typedef struct treenode_list_info treenode_list_info;
struct treenode_list_info{
  char             catalog_name[32];
  int              n_list;
  int              n_list_alloc;
  int              flag_groups_list;
  tree_node_info **list;
  ADaPS           *data;
};
// A --- Datatype definitions --- A

#ifdef __cplusplus
extern "C" {
#endif   
// V --- ANSI-C function definitions --- V
halo_properties_info *fetch_treenode_properties(tree_info *trees,tree_node_info *halo);
void   init_treenode_list(const char          *catalog_name,
                          int                  n_list_alloc,
                          treenode_list_info **list);
void   init_treenode_info_data(treenode_list_info  *list,
                               void               **rval,
                               SID_Datatype         data_type,
                               const char          *name,
                               ...);
void   free_treenode_list(treenode_list_info **list);
void   add_to_treenode_list(treenode_list_info *list,tree_node_info *node);
int    check_treenode_if_main_progenitor(tree_node_info *halo);
int    check_treenode_if_branch_start(tree_node_info *halo);
int    check_treenode_if_satellite(tree_node_info *halo);
int    check_treenode_if_central(tree_node_info *halo);
int    check_treenode_if_fragmented(tree_node_info *halo);
int    find_treenode_main_progenitor(tree_info *trees,tree_node_info *halo,tree_node_info **main_progenitor);
int    find_treenode_branch_root(tree_info *trees,tree_node_info *halo,tree_node_info **branch_root);
int    find_treenode_branch_leaf(tree_info *trees,tree_node_info *halo,tree_node_info **branch_leaf);
int    find_treenode_formation(tree_info       *trees,
                               tree_node_info  *halo,
                               double           f,
                               tree_node_info **peak_mass,
                               tree_node_info **fraction_of_peak_mass);
int    find_treenode_accretion(tree_info       *trees,
                               tree_node_info  *halo,
                               tree_node_info **first_satellite,
                               tree_node_info **join_current_group);
int    find_treenode_markers(tree_info *trees,tree_node_info *halo,tree_markers_info *markers);
int    fetch_treenode_snapshot(tree_info *trees,tree_node_info *halo);
int    fetch_treenode_snap_tree(tree_info *trees,tree_node_info *halo);
int    fetch_treenode_file_index(tree_info *trees,tree_node_info *halo);
double fetch_treenode_Mvir(tree_info *trees,tree_node_info *halo);
int    fetch_treenode_n_particles(tree_info *trees,tree_node_info *halo);
int    find_treesnap_z(tree_info *trees,double z_exact);
void   write_treenode_list_markers(tree_info *trees,treenode_list_info *list);
void   write_treenode_list_markers_header(tree_info *trees,treenode_list_info *list,FILE *fp_props_out);
void   write_tree_branches(tree_info *trees,tree_node_info **list_in,int n_list_in,int mode,const char *catalog_name);
void   write_tree_branch_ascii(tree_info *trees,tree_node_info *halo,const char *filename_out,const char *trees_name);
void   average_tree_branches(const char *catalog_name);
void   analyze_halos_and_N_subhalos(tree_info  *trees,
                                    const char *catalog_root,
                                    double      z_obs,
                                    double      M_cut_lo,
                                    double      M_cut_hi,
                                    int         n_subgroups_track_max);
// A --- ANSI-C function definitions --- A
#ifdef __cplusplus
}
#endif

#endif
