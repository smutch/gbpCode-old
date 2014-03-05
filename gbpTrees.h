#ifndef GBPTREES_AWAKE
#define GBPTREES_AWAKE
#include <gbpHalos.h>
#include <gbpCosmo.h>

#define K_MATCH_SUBGROUPS 0
#define K_MATCH_GROUPS    1

// Tree finalization modes
#define TREE_PROGENITOR_ORDER_N_PARTICLES 4
#define TREE_PROGENITOR_ORDER_DELUCIA     2
#define TREE_PROGENITOR_ORDER_DEFAULT     TREE_PROGENITOR_ORDER_N_PARTICLES
#define TREE_SUBSTRUCTURE_ORDER_DEFAULT   1
#define TREE_MODE_DEFAULT                 (TREE_SUBSTRUCTURE_ORDER_DEFAULT|TREE_PROGENITOR_ORDER_DEFAULT)

// If any of these are changed, don't forget to modify parse_match_type.c
#define TREE_CASE_SIMPLE                        1       // Set when a halo has a file_offset=1 
#define TREE_CASE_MAIN_PROGENITOR               2       // Set for the progenitor with the highest match score. (propagated for ghosts)
#define TREE_CASE_MERGER                        4       // Set when new IDs are created (ie. last point the halo was seen).
                                                        //    Set only for the last ghost in ghost-populated trees for mergers w/ offset>1.
#define TREE_CASE_DROPPED                       8       // Set if file_offset>1 and TREE_CASE_MATCHED_TO_BRIDGE is not set
#define TREE_CASE_STRAYED                       16      // Set for halos for which a descendant was not found
#define TREE_CASE_BRIDGED                       32      // Set for halos with multiple back-matches from halos with unique IDs
#define TREE_CASE_EMERGED_CANDIDATE             64      // Set when a halo is identified as a unique back-match to a halo marked TREE_CASE_BRIDGED 
                                                        //    and is not identified as the BRIDGE's main descendant
#define TREE_CASE_EMERGED                       128     // Set when a match is made identifying this halo as emerged
#define TREE_CASE_NO_PROGENITORS                256     // Set for halos that have no progenitors.
#define TREE_CASE_FRAGMENTED_STRAYED            512     // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                        //    decendant_id!=a valid id (ie they are not a progenitor of anything). (propagated for ghosts)
#define TREE_CASE_FRAGMENTED_RETURNED           1024    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                        //    decendant_id==the id of the halo they are emerged from. (propagated for ghosts)
#define TREE_CASE_FRAGMENTED_EXCHANGED          2048    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                        //    decendant_id!=the id of the halo they are emerged but is nevertheless valid 
                                                        //    (ie. they are still a progenitor of something). (propagated for ghosts)
#define TREE_CASE_MATCHED_TO_BRIDGE             4096    // Set when a halo is matched to one with TREE_CASE_BRIDGED set
#define TREE_CASE_BRIDGE_DEFAULT                8192    // Set when a halo matched to a bridge is not matched to any emerged candidate halos
#define TREE_CASE_GHOST                         16384   // Marks ghost halos in ghost-populated trees
#define TREE_CASE_GHOST_NULL                    32768   // Marks a ghost halo where a subgroup is it's own group.
                                                        //    This is a default behaviour that occurs when a group is strayed but one of 
                                                        //    it's subgroups isn't.
#define TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED 65536   // For internal use.  This should never be seen in the output.
#define TREE_CASE_BRIDGE_FINALIZE               131072  // For internal use.  This should never be seen in the output.
#define TREE_CASE_UNPROCESSED                   262144  // For internal use.  This should never be seen in the output.
#define TREE_CASE_INVALID                       524288  // For internal use.  This should never be seen in the output.
#define TREE_CASE_FRAGMENTED_NEW                (TREE_CASE_EMERGED_CANDIDATE+TREE_CASE_NO_PROGENITORS)

#define TREE_HORIZONTAL_READ_DEFAULT   0
#define TREE_HORIZONTAL_READ_EXTENDED  1
#define TREE_HORIZONTAL_STORE_EXTENDED 2
#define TREE_HORIZONTAL_STORE_GHOSTS   4

#define TREE_HORIZONTAL_WRITE_DEFAULT           0
#define TREE_HORIZONTAL_WRITE_ALLCASES          1
#define TREE_HORIZONTAL_WRITE_NOCASES           2
#define TREE_HORIZONTAL_WRITE_EXTENDED          4
#define TREE_HORIZONTAL_WRITE_GHOSTS            8
#define TREE_HORIZONTAL_WRITE_CHECK_FRAGMENTED 16

#define WRITE_MATCHES_MODE_TREES     0
#define WRITE_MATCHES_MODE_SINGLE    1
#define WRITE_MATCHES_CHECK_HEADER   2
#define WRITE_MATCHES_PERFORM_CHECK  4
#define WRITE_MATCHES_MODE_DEFAULT   WRITE_MATCHES_MODE_SINGLE

#define INIT_TREE_DATA_SUBGROUPS  0
#define INIT_TREE_DATA_GROUPS     1
#define INIT_TREE_DATA_DEFAULT    INIT_TREE_DATA_SUBGROUPS

#define READ_TREES_CATALOGS_SUBGROUPS   0
#define READ_TREES_CATALOGS_GROUPS      1
#define READ_TREES_CATALOGS_PROFILES    2
#define READ_TREES_CATALOGS_SHORT       4
#define READ_TREES_CATALOGS_ALL         READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS|READ_TREES_CATALOGS_PROFILES
#define READ_TREES_CATALOGS_DEFAULT     READ_TREES_CATALOGS_SUBGROUPS

#define READ_TREES_MATCH_SCORES_SUBGROUPS   0
#define READ_TREES_MATCH_SCORES_GROUPS      1
#define READ_TREES_MATCH_SCORES_ALL         READ_TREES_MATCH_SCORES_GROUPS|READ_TREES_MATCH_SCORES_SUBGROUPS
#define READ_TREES_MATCH_SCORES_DEFAULT     READ_TREES_MATCH_SCORES_SUBGROUPS

// Data structures for horizontal tree construction
typedef struct tree_horizontal_stats_info tree_horizontal_stats_info;
struct tree_horizontal_stats_info {
   int n_halos;
   int n_simple;
   int n_mergers;
   int n_strayed;
   int n_dropped;
   int n_bridged;
   int n_bridge_progenitors;
   int n_emerged;
   int n_fragmented_strayed;
   int n_fragmented_returned;
   int n_fragmented_exchanged;
   int n_emerged_progenitors;
   int n_invalid;
   int n_unprocessed;
   int max_strayed_size;
   int max_dropped_size;
   int max_bridged_size;
   int max_bridge_progenitor_size;
   int max_emerged_size;
   int max_fragmented_strayed_size;
   int max_fragmented_returned_size;
   int max_fragmented_exchanged_size;
   int max_emerged_progenitor_size;
   int max_id;
};

typedef struct tree_horizontal_info tree_horizontal_info;
typedef struct match_info bridge_info;
typedef struct match_info match_info;
struct match_info{
  tree_horizontal_info *halo;
  float                 score;
};
struct tree_horizontal_info{
  int          id;                 // This halo's id
  int          main_progenitor_id; // This halo's main progenitor id
  int          tree_id;            // This halo's tree id
  int          type;               // A bit-wise switch characterising this halo's matching
  int          file;               // This halo's snapshot index (ie. 0->n_snaps_used_in_trees-1)
  int          snap;               // This halo's snapshot number
  int          n_particles;        // Number of particles in this halo
  int          n_particles_parent; // Number of particles in this halo's parent halo
  int          n_bridges;          // The number of cnadidate emerged halos back-matched to this halo
  int          n_progenitors;      // The number of progenitors pointing to this halo
  int          index;              // This halo's index in the halo catalog
  match_info   first_progenitor;   // Pointer to this halo's first progenitor
  match_info   last_progenitor;    // Pointer to this halo's last  progenitor
  match_info   next_progenitor;    // Pointer to this halo's next  progenitor
  bridge_info *bridges;            // Contains the pointer information for all of the back-matches to this halo
  match_info   bridge_forematch;   // Pointer to a possible initial match to a bridged halo
  match_info   bridge_backmatch;   // Pointer to a possible back-matched bridged halo
  match_info   descendant;         // Contains all the needed pointers to the descendant
};

typedef struct tree_horizontal_extended_info tree_horizontal_extended_info;
struct tree_horizontal_extended_info{
  int          id;                 // This halo's id
  int          n_particles;        // Number of particles in this halo
  int          tree_id;            // This halo's tree id
  int          descendant_id;      // This halo's main progenitor id
  int          type;               // A bit-wise switch characterising this halo's matching
  int          file_offset;        // This halo's snapshot index (ie. 0->n_snaps_used_in_trees-1)
  int          n_particles_parent; // Number of particles in this halo's parent
  int          n_particles_desc;   // Number of particles in this halo's descendant
  int          n_particles_proj;   // Number of particles in this halo's progenitor
  int          score_desc;         // Matching score of this halo to it's descendant
  int          score_prog;         // Matching score of this halo to it's descendant
  int          snap_bridge;        // File index of any halo that this halo may be matched to
  int          index_bridge;       // Index of any bridge this halo may be matched to
  int          id_bridge;          // ID of any bridge this halo may be matched to
  int          index;              // This halo's index in the halo catalog
};

typedef struct tree_horizontal_ghost_interp_info tree_horizontal_ghost_interp_info;
struct tree_horizontal_ghost_interp_info{
  int    file_start;
  int    index_start;
  int    file_stop;
  int    index_stop;
  double time_start;
  double time_stop;
};

typedef struct tree_horizontal_ghost_subgroup_info tree_horizontal_ghost_subgroup_info;
struct tree_horizontal_ghost_subgroup_info{
  int          halo_index;                                // This halo's file index (populated after all the ghosts are added)
  int          id;                                        // This halo's id
  int          type;                                      // A bit-wise switch characterising this halo's matching
  int          descendant_id;                             // This halo's main progenitor id
  int          tree_id;                                   // This halo's tree id
  int          file_offset;                               // This halo's descendant file offset (pre-ghost-populating)
  int          file_index;                                // This halo's descendant index
  tree_horizontal_ghost_interp_info    interp;
  tree_horizontal_ghost_subgroup_info *descendant;        // Used to figure out 'file_index' once ghosts have been created for the next snap
  tree_horizontal_ghost_subgroup_info *next_substructure; // Used for substructure linked lists
};

typedef struct tree_horizontal_ghost_group_info tree_horizontal_ghost_group_info;
struct tree_horizontal_ghost_group_info{
  int          halo_index;                                 // This halo's file index (populated after all the ghosts are added)
  int          id;                                         // This halo's id
  int          type;                                       // A bit-wise switch characterising this halo's matching
  int          descendant_id;                              // This halo's main progenitor id
  int          tree_id;                                    // This halo's tree id
  int          file_offset;                                // This halo's descendant file offset (pre-ghost-populating)
  int          file_index;                                 // This halo's descendant index
  int          n_subgroups;                                // Number of substructures in this group
  tree_horizontal_ghost_interp_info    interp;
  tree_horizontal_ghost_subgroup_info *first_substructure; // Used for substructure linked lists
  tree_horizontal_ghost_subgroup_info *last_substructure;  // Used for substructure linked lists
};

// Data structures for vertical tree construction
typedef struct tree_vertical_node_info tree_vertical_node_info;
struct tree_vertical_node_info{
  halo_info                halo;
  int                      depth_first_index;
  int                      group_id;
  int                      halo_id;
  int                      descendant_id;
  int                      n_progenitors;
  tree_vertical_node_info *descendant;
  tree_vertical_node_info *progenitor_first;
  tree_vertical_node_info *progenitor_next;
  tree_vertical_node_info *progenitor_last;
  tree_vertical_node_info *group_halo_first;
  tree_vertical_node_info *group_halo_next;
  tree_vertical_node_info *neighbour_halo_next;
  tree_vertical_node_info *next; // Points to the next halo, in the order they are added to the tree
};

typedef struct tree_vertical_info tree_vertical_info;
struct tree_vertical_info{
  tree_vertical_node_info  *root;
  tree_vertical_node_info  *last_leaf;
  int                      *n_neighbours;
  tree_vertical_node_info **neighbour_halos;
  tree_vertical_node_info **neighbour_halo_last;
};

// Data structures for tree analysis
typedef struct tree_node_info tree_node_info;
struct tree_node_info{
  int             n_progenitors;
  int             n_substructures;
  int             snap_tree;
  int             file_index;
  int             neighbour_index;
  int             n_particles;
  float           match_score;
  int             tree_case;
  // Pointers for the substructure heirarchy
  tree_node_info *parent;
  tree_node_info *substructure_first;   // for substructure in this halo's parent
  tree_node_info *substructure_last;    // for substructure in this halo's parent
  tree_node_info *substructure_next;    // for substructure in this halo's parent
  // Merger tree pointers
  tree_node_info *descendant;
  tree_node_info *progenitor_first;
  tree_node_info *progenitor_last;
  tree_node_info *progenitor_next;
  // Bulk processing pointers
  tree_node_info *next_neighbour;  // This halo's snapshot
  tree_node_info *next_in_forest;  // This halo's forest
};

typedef struct tree_info tree_info;
struct tree_info{
  // Counts etc
  int              i_read_start;
  int              i_read_stop;
  int              i_read_step;
  int              n_search;
  int              n_snaps;
  int              n_forests;
  int              n_forests_local;
  int             *snap_list;
  double          *z_list;
  double          *t_list;
  int             *n_groups_snap_local;
  int             *n_subgroups_snap_local;
  int             *n_groups_forest_local;
  int             *n_subgroups_forest_local;
  int              max_n_groups_snap_local;
  int              max_n_subgroups_snap_local;
  int              max_n_groups_forest_local;
  int              max_n_subgroups_forest_local;
  // Pointers
  tree_node_info **first_neighbour_groups;
  tree_node_info **first_neighbour_subgroups;
  tree_node_info **last_neighbour_groups;
  tree_node_info **last_neighbour_subgroups;
  tree_node_info **first_in_forest_groups;
  tree_node_info **first_in_forest_subgroups;
  tree_node_info **last_in_forest_groups;
  tree_node_info **last_in_forest_subgroups;
  // A place to add other data to
  ADaPS           *data;
};

// Structure needed for passing parameters to the tree data ADaPS deallocation function
typedef struct store_tree_data_free_parms_info store_tree_data_free_parms_info;
struct store_tree_data_free_parms_info{
  int n_snaps;
};

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
void read_trees(char       *filename_tree_root,
                char       *filename_cat_root,
                int         mode_progenitor,
                tree_info **trees);
void init_tree_data(tree_info    *trees,
                    void       ***rval,
                    size_t        data_size,
                    int           mode,
                    const char   *name,
                    ...);
void free_tree_data(void **tree_data,void *params);
void read_trees_match_scores(tree_info *trees,
                             char      *filename_SSimPL_dir,
                             int        mode);
void read_trees_catalogs(tree_info              *trees,
                         char                   *filename_SSimPL_dir,
                         char                   *filename_catalog_name,
                         int                     mode);
void read_AHF_for_trees(char       *filename_root,
                        int         i_file,
                        plist_info *plist,
                        char       *catalog_name,
                        int         mode);
void read_matches(char    *filename_root_matches,
                  int      i_read,
                  int      j_read,
                  int      n_halos_max,
                  int      mode,
                  int     *n_groups_i,
                  int     *n_groups_j,
                  int     *n_particles_i_in,
                  int     *n_particles_j_in,
                  int     *n_sub_group_i,
                  int     *n_sub_group_j,
                  int     *match_ids,
                  float   *match_score,
                  size_t  *match_index);
int check_for_matching_input_files(const char *filename_root_in,int i_read);
int check_goodness_of_match(int n_particles_in,float match_score);
int compute_single_matches(char   *filename_root_in,
                           char   *filename_root_out,
                           int     i_read_1,
                           int     i_read_2);
int compute_trees_matches(char   *filename_halo_root_in,
                          char   *filename_root_out,
                          int     i_read_stop,
                          int     i_read_start,
                          int     i_read_step,
                          int     n_search,
                          int     mode);
void write_match_results(char       *filename_out_dir,
                         char       *filename_out_root,
                         int         i_read,
                         int         j_read,
                         plist_info *plist1,
                         plist_info *plist2,
                         int         k_match,
                         int         mode);
void compute_trees_horizontal_stats(void *halos_in,int n_halos,int n_halos_max,tree_horizontal_stats_info *stats,int flag_write_cases);

void compute_trees_horizontal_fragmented(int         *n_groups,
                                         int         *n_subgroups,
                                         int        **n_subgroups_group,
                                         int          i_file_start,
                                         int          i_write_last,
                                         int          j_write_last,
                                         int          l_write_last,
                                         int          i_read_stop,
                                         int          i_read_step,
                                         int          max_tree_id_subgroup,
                                         int          max_tree_id_group,
                                         int          n_subgroups_max,
                                         int          n_groups_max,
                                         int          n_search,
                                         int          n_files,
                                         int          n_wrap,
                                         int          n_k_match,
                                         double      *a_list,
                                         cosmo_info **cosmo,
                                         char        *filename_output_dir);
void compute_trees_horizontal_ghosts(int         *n_groups,
                                     int         *n_subgroups,
                                     int        **n_subgroups_group,
                                     int          i_read_start,
                                     int          i_file_start,
                                     int          i_write_last,
                                     int          j_write_last,
                                     int          l_write_last,
                                     int          i_read_stop,
                                     int          i_read_step,
                                     int          max_tree_id_subgroup,
                                     int          max_tree_id_group,
                                     int          n_subgroups_max,
                                     int          n_groups_max,
                                     int          n_search,
                                     int          n_files,
                                     int          n_wrap,
                                     int          n_k_match,
                                     double      *a_list,
                                     cosmo_info **cosmo,
                                     char        *filename_cat_root_in,
                                     char        *filename_output_dir);

int   set_match_id(match_info *match);
int   set_match_file(match_info *match);
int   set_match_snapshot(match_info *match);
int   set_match_type(match_info *match);
int   set_match_n_particles(match_info *match);
int   set_match_index(match_info *match);
float set_match_score(match_info *match);
void check_for_fragmented_halos(int k_match,tree_horizontal_info **groups,int n_groups,
                                int i_write,int j_write,int l_write,int n_wrap);
void add_to_trees_horizontal_stats(tree_horizontal_stats_info *stats,int id,int type,int n_particles);
void init_trees_horizontal_roots(tree_horizontal_info **groups,
                                 tree_horizontal_info **subgroups,
                                 int    *match_id,
                                 float  *match_score,
                                 size_t *match_index,
                                 int    *n_particles_groups,
                                 int    *n_particles_subgroups,
                                 int   **n_subgroups_group,
                                 int     n_groups_max,
                                 int     n_subgroups_max,
                                 char   *filename_root_matches,
                                 int     i_read_start,
                                 int     i_read_stop,
                                 int     i_read_step,
                                 int     i_file_start,
                                 int     n_wrap,
                                 int     n_halos_max,
                                 int    *max_id_group,
                                 int    *max_tree_id_group,
                                 int    *max_id_subgroup,
                                 int    *max_tree_id_subgroup);
void init_trees_horizontal_snapshot(tree_horizontal_info *halos_i,int n_halos_i,int i_read,int i_file,int n_halos_max);
void identify_bridges(tree_horizontal_info **halos,
                      tree_horizontal_info  *halos_i,
                      int     n_halos_i,
                      int    *match_id,
                      float  *match_score,
                      size_t *match_index,
                      int    *n_particles,
                      int     i_file,
                      int     i_read,
                      int     i_read_start,
                      int     i_read_stop,
                      int     i_read_step,
                      int     n_search,
                      int     n_wrap,
                      int     n_halos_max,
                      int     n_files,
                      char   *filename_root_matches,
                      int     flag_match_subgroups);
void construct_progenitors(tree_horizontal_info **halos,
                           tree_horizontal_info  *halos_i,
                           int   **n_subgroups_group,
                           int     n_halos_i,
                           int    *match_id,
                           float  *match_score,
                           size_t *match_index,
                           int    *n_particles,
                           int     i_file,
                           int     i_read,
                           int     i_read_start,
                           int     i_read_stop,
                           int     i_read_step,
                           int     n_search,
                           int     n_wrap,
                           int     n_halos_max,
                           int     n_files,
                           int     flag_fix_bridges,
                           int    *max_id,
                           int    *n_halos_1_matches,
                           int    *n_halos_2_matches,
                           char   *filename_root_matches,
                           char   *group_text_prefix,
                           int     flag_match_subgroups);
void clean_emerged_halo_list(tree_horizontal_info *halos_i,
                             int                   n_halos_i,
                             int                   i_file,
                             int                   n_search,
                             int                   n_files);
void apply_horizontal_tree_defaults(int                    n_halos_1_matches,
                                    int                    n_halos_i,
                                    tree_horizontal_info **halos,
                                    tree_horizontal_info  *halos_i,
                                    int                    i_file,
                                    int                    n_wrap,
                                    int                   *max_id,
                                    int                   *max_tree_id);
void init_trees_horizontal_stats(tree_horizontal_stats_info *stats,int n_halos);
void change_horizontal_ID_recursive(tree_horizontal_info *halo,int id_1,int id_2);
void write_trees_horizontal_log_file(char *filename_log,
                                     int   l_write,
                                     int   j_write,
                                     int   i_k_match,
                                     int   n_k_match,
                                     tree_horizontal_stats_info  *stats,
                                     double                      *a_list,
                                     cosmo_info                 **cosmo,
                                     int                          flag_init);
void write_trees_horizontal_report(int                   n_halos_i,
                                   int                   n_halos_max,
                                   tree_horizontal_info *halos_i);
void write_trees_horizontal(void **groups_in,   
                            void **subgroups_in,
                            int    n_groups,    int n_groups_max,   
                            int    n_subgroups, int n_subgroups_max,
                            int   **n_subgroups_group,
                            int     max_tree_id_subgroup,
                            int     max_tree_id_group,
                            int     i_write,
                            int     j_write,
                            int     l_write,
                            int     n_step,
                            int     n_search,
                            int     n_wrap,
                            int     i_file_start,
                            char   *filename_cat_root_in,
                            char   *filename_output_dir,
                            double *a_list,
                            cosmo_info **cosmo,
                            int     n_k_match,
                            int     flag_init_write,
                            int     mode);
void write_trees_horizontal_emerged_candidates(int                   i_read,
                                               int                   n_halos_i,
                                               tree_horizontal_info *halos_i,
                                               char                 *group_text_prefix,
                                               char                 *filename_output_dir,
                                               int                   flag_start_new_file);
void read_forests(char  *filename_root_in,
                  int    n_trees_group,
                  int    n_trees_subgroup,
                  int   *n_forests_group,
                  int   *n_forests_subgroup,
                  int   *n_forests_group_local,
                  int   *n_forests_subgroup_local,
                  int  **i_forest_group,
                  int  **i_forest_subgroup,
                  int  **n_halos_forest_group,
                  int  **n_halos_forest_subgroup,
                  int   *n_trees_forest_groups_max,
                  int   *n_trees_forest_subgroups_max,
                  int   *forest_lo_group_local,
                  int   *forest_hi_group_local,
                  int   *forest_lo_subgroup_local,
                  int   *forest_hi_subgroup_local,
                  int   *n_groups_local,
                  int   *n_subgroups_local,
                  int   *n_groups_max_snap_local,
                  int   *n_subgroups_max_snap_local);
void read_tree_final_totals(char *filename_output_dir_horizontal_trees,
                            int   i_read_start,
                            int   i_read_stop, 
                            int   i_read_step,
                            int  *i_read_last,
                            int  *n_snap,
                            int  *n_groups_max_in,
                            int  *n_subgroups_max_in,
                            int  *n_progenitors_max,
                            int  *n_trees_subgroup,
                            int  *n_trees_group);
void split_forests_n_ways(int  *n_halos_forest,
                          int   n_forests,
                          int   n_split,
                          int **tree_count_split,
                          int **forest_lo_split,
                          int **forest_hi_split);
int read_matches_header(char   *filename_root_in,
                        int     i_read_start,
                        int     i_read_stop,
                        int     i_read_step,
                        int    *n_files_return,
                        int   **n_subgroups_return,
                        int   **n_groups_return,
                        int    *n_subgroups_max,
                        int    *n_groups_max,
                        int    *n_halos_max);
void read_trees_horizontal(void **groups,   int *n_groups_in,
                           void **subgroups,int *n_subgroups_in,
                           int   *n_subgroups_group,
                           int   *n_trees_subgroup_in,
                           int   *n_trees_group_in,
                           int    i_read, // tree snapshot index
                           int    j_read, // actual snapshot index
                           int    l_read,
                           int    n_wrap,
                           char  *filename_output_dir,
                           int    mode);
void propagate_fragmented_halos(tree_horizontal_extended_info **groups,   int *n_groups,
                                tree_horizontal_extended_info **subgroups,int *n_subgroups,
                                int        **n_subgroups_group,
                                int          i_read, // tree snapshot index
                                int          j_read,
                                int          l_read,
                                int          i_read_step,
                                int          n_wrap);
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_wrap);
void process_ghosts(tree_horizontal_ghost_group_info    **groups,
                    tree_horizontal_ghost_subgroup_info **subgroups,
                    int        *n_groups,
                    int        *n_subgroups,
                    int       **n_subgroups_group,
                    int        *n_group_ghosts,
                    int        *n_subgroup_ghosts,
                    int        *n_group_ghosts_used,
                    int        *n_subgroup_ghosts_used,
                    int         i_read,
                    int         l_read,
                    int         j_file,
                    int         i_file_start,
                    int         n_search,
                    int         n_wrap,
                    int         n_files,
                    double     *a_list,
                    cosmo_info **cosmo);
void add_substructure_to_horizontal_tree_group(tree_horizontal_ghost_group_info    *group,
                                               tree_horizontal_ghost_subgroup_info *subgroup_descendant,
                                               tree_horizontal_ghost_subgroup_info *subgroup);
void read_catalog_ghost_interpolation(tree_horizontal_ghost_group_info     *groups,
                                      halo_properties_info                **group_properties,
                                      int                                   n_groups,
                                      tree_horizontal_ghost_subgroup_info  *subgroups,
                                      halo_properties_info                **subgroup_properties,
                                      int                                   n_subgroups,
                                      char                                 *filename_cat_root_in,
                                      int                                   i_read,
                                      int                                   j_read,
                                      int                                   n_wrap);
void write_ghost_catalog(tree_horizontal_ghost_group_info      *groups,
                         halo_properties_info                ***group_properties,
                         int                                    n_group_ghosts,
                         int                                    n_groups,
                         tree_horizontal_ghost_subgroup_info   *subgroups,
                         halo_properties_info                ***subgroup_properties,
                         int                                    n_subgroup_ghosts,
                         int                                    n_subgroups,
                         char                                  *filename_output_dir,
                         char                                  *filename_cat,
                         int                                    i_read,
                         int                                    j_read,
                         int                                    l_read,
                         int                                    n_wrap,
                         double                                *a_list,
                         cosmo_info                           **cosmo);
void compute_trees_horizontal(char   *filename_halos_root_in,
                              char   *filename_cat_root_in,
                              char   *filename_root_matches,
                              char   *filename_root_out,
                              double *a_list,
                              cosmo_info **cosmo,
                              int     i_read_start,
                              int     i_read_stop,
                              int     i_read_step,
                              int     n_search,
                              int     flag_fix_bridges,
                              int    *flag_clean);
void compute_forests(char *filename_root_out,int n_search_forests);
void compute_trees_vertical_old(char *filename_root_out,
                                char *filename_cat_root_in,
                                char *filename_snap_list_in,
                                int   n_files_groups,
                                int   n_files_subgroups,
                                int   n_search_forests,
                                int  *flag_clean);
void compute_trees_vertical(char *filename_SSimPL_dir,
                            char *filename_halo_version_root,
                            char *filename_trees_name,
                            int   n_dim_files);
void finalize_trees(tree_info *trees,
                    int        progenitor_mode);
void finalize_trees_vertical(tree_vertical_info **trees,
                             int                 *n_halos_tree,
                             int                  n_trees,
                             int                  n_snaps,
                             int                  progenitor_mode);
void read_tree_run_parameters(char *filename_root_out,
                              int  *i_read_start,
                              int  *i_read_stop,
                              int  *i_read_step,
                              int  *n_search,
                              int  *flag_fix_bridges,
                              int  *flag_compute_fragmented,
                              int  *flag_compute_ghosts);
void write_tree_run_parameters(char *filename_root_out,
                               int   i_read_start,
                               int   i_read_stop,
                               int   i_read_step,
                               int   n_search,
                               int   flag_fix_bridges,
                               int   flag_compute_fragmented,
                               int   flag_compute_ghosts);
void write_trees_vertical(tree_vertical_info **trees,
                          int                 *n_halos_tree_local,
                          int                  n_trees_local,
                          int                 *tree_lo_file,
                          int                 *tree_hi_file,
                          int                 *n_halos_file,
                          int                  n_files,
			  const char          *filename_root_out,
                          const char          *group_text_prefix);
void write_a_list(const char *filename_snap_list_in,
                  const char *filename_root_out,
                  int         i_read_start,
                  int         i_read_stop,
                  int         i_read_step);
void compute_trees_auxiliary(char *filename_root,
                             char *filename_snapshot_root,
                             char *filename_root_out,
                             int   i_read_start,
                             int   i_read_stop,
                             int   i_read_step,
                             int   n_search,
                             int   n_files_groups,
                             int   n_files_subgroups,
                             int  *flag_clean);
void init_trees(int         i_read_start,
                int         i_read_stop,
                int         i_read_step,
                int         n_search,
                int         n_forests,
                int         n_forests_local,
                tree_info **tree);
void free_trees(tree_info **tree);
int add_node_to_trees(tree_info        *trees,
                      int               i_forest,
                      int               tree_case,
                      int               n_particles,
                      int               halo_snap,
                      int               halo_index,
                      int               descendant_snap,
                      int               descendant_index,
                      int             **halo_indices,
                      tree_node_info ***halo_array,
                      int               n_wrap,
                      tree_node_info   *group_node,
                      tree_node_info  **new_node);
void init_trees_vertical(int n_snaps,tree_vertical_info **tree);
void free_trees_vertical(tree_vertical_info **tree);
int  add_node_to_vertical_tree(tree_vertical_info  *tree,
                               int                  match_type,
                               int                  halo_id,
                               int                  group_id,
                               int                  descendant_id,
                               int                  halo_snap,
                               int                  descendant_snap,
                               halo_info           *properties);

void assign_group_order(tree_info *tree,int mode);
void assign_progenitor_order(tree_info *tree,int mode);
void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode);
void assign_depth_first_index_recursive(tree_node_info *tree,int *depth_first_index);
void assign_unique_tree_ids_recursive(tree_node_info *tree_node,int i_tree);
int  construct_unique_tree_id(tree_node_info *tree_node,int tree_number);
void compute_progenitor_score_recursive(tree_node_info *tree,int *M_i,int mode);
void compute_substructure_order_recursive(tree_node_info *parent,int *score_parent,int mode);
void compute_progenitor_order_recursive(tree_node_info *descendant,int *score_descendant,int mode);
void compute_trees_analysis(tree_info *trees,char *filename_out_root);

void assign_group_subgroup_order_vertical(tree_vertical_info *tree,int i_snap,int mode);
void assign_progenitor_order_vertical_recursive(tree_vertical_node_info *tree,int *M_i,int mode);
void assign_depth_first_index_vertical_recursive(tree_vertical_node_info *tree,int *depth_first_index);
void assign_unique_vertical_tree_ids_recursive(tree_vertical_node_info *tree_node,int i_tree);
int  construct_unique_vertical_tree_id(tree_vertical_node_info *tree_node,int tree_number);
void compute_progenitor_score_vertical_recursive(tree_vertical_node_info *tree,int *M_i,int mode);

#ifdef __cplusplus
}
#endif

#endif

