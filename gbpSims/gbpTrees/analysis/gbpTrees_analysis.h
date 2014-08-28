#ifndef GBPTREES_ANALYSIS_AWAKE
#define GBPTREES_ANALYSIS_AWAKE

// V Preprocessor definitions V
#define WRITE_TREENODE_LIST_PROPERTIES_N 8

#define COMPUTE_ACCRETION_ANALYSIS_GROUPS    TTTP01
#define COMPUTE_ACCRETION_ANALYSIS_SUBGROUPS TTTP02
#define COMPUTE_ACCRETION_ANALYSIS_DEFAULT   COMPUTE_ACCRETION_ANALYSIS_GROUPS
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
  tree_node_info *merger_33pc_remnant;
  tree_node_info *merger_33pc_host;
  tree_node_info *merger_33pc_merger;
  tree_node_info *merger_10pc_remnant;
  tree_node_info *merger_10pc_host;
  tree_node_info *merger_10pc_merger;
  double          M_peak;
};

typedef struct tree_markers_stats_info tree_markers_stats_info;
struct tree_markers_stats_info{
  double t_half_peak_mass_ranges[2][2];
  double t_half_peak_mass_peak;
  double t_merger_33pc_ranges[2][2];
  double t_merger_33pc_peak;
  double t_merger_10pc_ranges[2][2];
  double t_merger_10pc_peak;
  double M_peak_ranges[2][2];
  double M_peak_peak;
  double M_vir_ranges[2][2];
  double M_vir_peak;
};

typedef struct treenode_list_info treenode_list_info;
struct treenode_list_info{
  char             catalog_name[MAX_FILENAME_LENGTH];
  int              n_list;
  int              n_list_local;
  int              n_list_alloc;
  int              flag_groups_list;
  tree_node_info **list;
  ADaPS           *data;
};

#define TREENODE_HIST_LOG_X       1
#define TREENODE_HIST_LOG_Y       2
#define TREENODE_HIST_DEFAULT     0
#define TREENODE_HIST_N_ARGS_MAX  3
#define TREENODE_HIST_N_PROPS     3
#define TREENODE_HIST_NAME_LENGTH 128
typedef struct treenode_hist_info treenode_hist_info;
struct treenode_hist_info{
  char    name[TREENODE_HIST_NAME_LENGTH];
  int     x_prop;
  int     y_prop;
  double  x_args_d[TREENODE_HIST_N_ARGS_MAX];
  double  y_args_d[TREENODE_HIST_N_ARGS_MAX];
  int     x_args_i[TREENODE_HIST_N_ARGS_MAX];
  int     y_args_i[TREENODE_HIST_N_ARGS_MAX];
  int     n_x;
  int     n_y;
  int     flag_log_x;
  int     flag_log_y;
  int    *array;
};
typedef struct treenode_hist_props_info treenode_hist_props_info;
struct treenode_hist_props_info{
   char         name[TREENODE_HIST_N_PROPS][TREENODE_HIST_NAME_LENGTH];
   char         axis_text[TREENODE_HIST_N_PROPS][TREENODE_HIST_NAME_LENGTH];
   char         log_axis_text[TREENODE_HIST_N_PROPS][TREENODE_HIST_NAME_LENGTH];
   int          n_args[TREENODE_HIST_N_PROPS];
   SID_Datatype arg_type[TREENODE_HIST_N_PROPS][TREENODE_HIST_N_ARGS_MAX];
};
#ifdef _MAIN
treenode_hist_props_info treenode_hist_props={{"z","M","N"},
                                              {"$z$","$M_{\\rm{vir}}[h^{-1} M_\\odot]$","$\\rm{n_p}$"},
                                              {"$\\log_{10}{z}$","$\\log_{10}{\\rm{M}_{\\rm{vir}}[h^{-1} \\rm{M}_\\odot]}$","$\\log_{10}{\\rm{n_p}}$"},
                                              {  1,  3,  3},
                                              {{SID_INT,   SID_INT,   SID_INT},
                                               {SID_DOUBLE,SID_DOUBLE,SID_INT},
                                               {SID_DOUBLE,SID_DOUBLE,SID_INT}}};
#else
extern treenode_hist_props_info treenode_hist_props;
#endif

// A --- Datatype definitions --- A

#ifdef __cplusplus
extern "C" {
#endif   
// V --- ANSI-C function definitions --- V
halo_properties_info *fetch_treenode_properties(tree_info *trees,tree_node_info *halo);
void   init_treenode_hist(tree_info           *trees,
                          const char          *hist_name,
                          const char          *x_name_in,
                          const char          *y_name_in,
                          int                  mode,
                          treenode_hist_info **hist, ...);
void   add_to_treenode_hist(tree_info           *trees,
                            treenode_hist_info  *hist,
                            tree_node_info      *current_halo);
void   free_treenode_hist  (treenode_hist_info **hist);
void   init_treenode_list(const char          *catalog_name,
                          int                  n_list_alloc,
                          treenode_list_info **list);
void   init_treenode_info_data(treenode_list_info  *list,
                               void               **rval,
                               SID_Datatype         data_type,
                               const char          *name,
                               ...);
void finalize_treenode_list(tree_info *trees,treenode_list_info *list);
void   free_treenode_list(treenode_list_info **list);
void   add_to_treenode_list(treenode_list_info *list,tree_node_info *node);
int    check_treenode_if_snap_equals_given(tree_node_info *halo,int snap_tree_given);
int    check_treenode_if_main_progenitor(tree_node_info *halo);
int    check_treenode_if_merger(tree_node_info *halo);
int    check_treenode_if_branch_start(tree_info *trees,tree_node_info *halo);
int    check_treenode_if_satellite(tree_node_info *halo);
int    check_treenode_if_central(tree_node_info *halo);
int    check_treenode_if_fragmented(tree_node_info *halo);
int    find_treenode_snap_equals_given(tree_info *trees,tree_node_info *halo,int snap_tree_given,tree_node_info **treenode_return);
int    find_treenode_main_progenitor(tree_info *trees,tree_node_info *halo,tree_node_info **main_progenitor);
int    find_treenode_branch_root(tree_info *trees,tree_node_info *halo,tree_node_info **branch_root);
int    find_treenode_branch_leaf(tree_info *trees,tree_node_info *halo,tree_node_info **branch_leaf);
int    find_treenode_formation(tree_info       *trees,
                               tree_node_info  *halo,
                               double           f,
                               tree_node_info **peak_mass,
                               tree_node_info **fraction_of_peak_mass);
int find_treenode_last_merger(tree_info       *trees,
                              tree_node_info  *halo,
                              double           fraction,
                              tree_node_info **remnant,
                              tree_node_info **merger_host,
                              tree_node_info **merged_halo);
int    find_treenode_accretion(tree_info       *trees,
                               tree_node_info  *halo,
                               tree_node_info **first_satellite,
                               tree_node_info **join_current_group);
int    find_treenode_markers(tree_info *trees,tree_node_info *halo,tree_markers_info *markers);
int    find_treenode_markers_all(tree_info *trees,tree_markers_info ***markers,int mode);
int    init_treenode_markers_all(tree_info *trees,tree_markers_info ***markers,int mode);
int    free_treenode_markers_all(tree_info *trees,tree_markers_info ***markers);
void   write_treenode_markers_all(tree_info *trees,tree_markers_info **markers,char *filename_output_root,int mode);
int    find_treenode_markers_recursive(tree_info          *trees,
                                       tree_markers_info **markers_array,
                                       tree_node_info     *halo,
                                       int                 flag_is_main_progenitor,
                                       double             *M_return,
                                       tree_markers_info **markers_descendant);
void compute_treenode_list_marker_stats(tree_info                *trees,
                                        treenode_list_info       *list,
                                        tree_markers_info       **markers_all,
                                        tree_markers_stats_info  *stats,
                                        int                     **n_hist_count,
                                        int                      *n_hist);
void   compute_trees_analysis                (tree_info *trees);
void   compute_trees_emerged_halo_analysis   (tree_info *trees,char *filename_out_root);
void   compute_trees_fragmented_halo_analysis(tree_info *trees,char *filename_out_root);
void   compute_trees_merger_analysis         (tree_info *trees,char *filename_out_root);
void   compute_trees_dropped_halo_analysis   (tree_info *trees,char *filename_out_root);
void   compute_trees_strayed_halo_analysis   (tree_info *trees,char *filename_out_root);
void   compute_marker_analysis(tree_info *trees,tree_markers_info **markers,const char *catalog_root,const char *filename_out_root,int mode);

double fetch_treenode_delta_t(tree_info *trees,tree_node_info *halo_1,tree_node_info *halo_2);
int    fetch_treenode_snapshot(tree_info *trees,tree_node_info *halo);
int    fetch_treenode_snap_tree(tree_info *trees,tree_node_info *halo);
int    fetch_treenode_file_index(tree_info *trees,tree_node_info *halo);
double fetch_treenode_Mvir(tree_info *trees,tree_node_info *halo);
double fetch_treenode_x(tree_info *trees,tree_node_info *halo);
double fetch_treenode_y(tree_info *trees,tree_node_info *halo);
double fetch_treenode_z(tree_info *trees,tree_node_info *halo);
double fetch_treenode_vx(tree_info *trees,tree_node_info *halo);
double fetch_treenode_vy(tree_info *trees,tree_node_info *halo);
double fetch_treenode_vz(tree_info *trees,tree_node_info *halo);
int    fetch_treenode_n_particles(tree_info *trees,tree_node_info *halo);
int    find_treesnap_z(tree_info *trees,double z_exact);
void   write_treenode_list_markers           (tree_info *trees,const char *filename_out_root,treenode_list_info *list);
void   write_treenode_list_markers_header    (tree_info *trees,treenode_list_info *list,FILE *fp_props_out);
void   write_treenode_list_data              (tree_info *trees,const char *filename_out_root,treenode_list_info *list);
void   write_treenode_list_data_header       (tree_info *trees,treenode_list_info *list,FILE *fp_props_out);
void   write_treenode_list_properties        (tree_info *trees,const char *filename_out_root,treenode_list_info *list);
void   write_treenode_list_properties_header (tree_info *trees,treenode_list_info *list,FILE *fp_props_out);
void   write_treenode_list_properties_set_ith(tree_info *trees,tree_node_info *current_halo,int i_write,char *data_name,SID_Datatype *data_type,int *data_i,double *data_d);
void   write_tree_branches                   (tree_info *trees,tree_node_info **list_in,int n_list_in,int mode,const char *filename_out_dir,const char *catalog_name);
void   write_tree_branch_ascii               (tree_info *trees,tree_node_info *halo,const char *filename_out,const char *trees_name);
void   write_treenode_all_hist               (tree_info *trees,const char *filename_root_out);
void   write_treenode_hist(tree_info          *trees,
                           const char         *filename_out_root,
                           treenode_hist_info *hist);
void   write_treenode_list_hist(tree_info          *trees,
                                const char         *filename_root_out,
                                treenode_list_info *list);
void   average_tree_branches(const char *catalog_name);
void   analyze_halos_and_N_subhalos(tree_info  *trees,
                                    const char *filename_out_root,
                                    const char *catalog_root,
                                    double      z_obs_exact,
                                    double      M_cut_lo,
                                    double      M_cut_hi,
                                    int         n_subgroups_track_max);
// A --- ANSI-C function definitions --- A
#ifdef __cplusplus
}
#endif

#endif
