#ifndef GBPTREES_AWAKE
#define GBPTREES_AWAKE
#include <gbpHalos.h>

#define TREE_PROGENITOR_ORDER_DEFAULT 0
#define TREE_PROGENITOR_ORDER_DELUCIA 2

#define TREE_CASE_SIMPLE                        1
#define TREE_CASE_UNPROCESSED                   2
#define TREE_CASE_INVALID                       4
#define TREE_CASE_STRAYED                       8
#define TREE_CASE_SPUTTERED                     16
#define TREE_CASE_DROPPED                       32
#define TREE_CASE_FOUND                         64
#define TREE_CASE_MERGER                        128
#define TREE_CASE_BRIDGED                       256
#define TREE_CASE_EMERGED                       512
#define TREE_CASE_BRIDGE_PROGENITOR             1024
#define TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED 4096
#define TREE_CASE_BRIDGE_FINALIZE               8192

typedef struct tree_horizontal_stats_info tree_horizontal_stats_info;
struct tree_horizontal_stats_info {
   int n_halos;
   int n_simple;
   int n_mergers;
   int n_strayed;
   int n_sputtered;
   int n_dropped;
   int n_bridged;
   int n_bridge_progenitors;
   int n_emerged;
   int n_emerged_lost;
   int n_emerged_progenitors;
   int max_strayed_size;
   int max_sputtered_size;
   int max_dropped_size;
   int max_bridged_size;
   int max_bridge_progenitor_size;
   int max_emerged_size;
   int max_emerged_lost_size;
   int max_emerged_found_diff;
   int max_emerged_found_diff_size;
   int max_emerged_progenitor_size;
   int max_id;
};

typedef struct match_info match_info;
struct match_info{
  int    id;
  int    file;
  float  score;
  size_t index;
  int    n_particles;
};

typedef struct tree_horizontal_info tree_horizontal_info;
struct tree_horizontal_info{
  int         id;            // This halo's id
  int         tree_id;       // This halo's tree id
  match_info  descendant;    // Contains all the needed pointers to the descendant
  int         n_progenitors; // The number of progenitors pointing to this halo
  int         type;          // A flag list characterising this halo's matching
  int         n_particles;   // Number of particles in this halo
  int         n_bridges;     // The number of bridges back-matched to this halo
  match_info *bridges;       // Contains the pointer information for all of the back-matches to this halo
  match_info *bridge_match;  // In the event that this halo gets matched to a bridge initially,
                             //   this holds the information about the bridge for statistics writing
};

typedef struct tree_node_info tree_node_info;
struct tree_node_info{
  halo_info       halo;
  int             depth_first_index;
  int             group_id;
  int             halo_id;
  int             n_progenitors;
  tree_node_info *descendant;
  tree_node_info *progenitor_first;
  tree_node_info *progenitor_next;
  tree_node_info *progenitor_last;
  tree_node_info *group_halo_first;
  tree_node_info *group_halo_next;
  tree_node_info *neighbour_halo_next;
  tree_node_info *next; // Points to the next halo, in the order they are added to the tree
};

typedef struct tree_info tree_info;
struct tree_info{
  tree_node_info  *root;
  tree_node_info  *last_leaf;
  int             *n_neighbours;
  tree_node_info **neighbour_halos;
  tree_node_info **neighbour_halo_last;
};

void read_trees(char             *filename_trees_root,
                int               i_file_start,
                int               i_file_stop,
                int               mode,
                tree_node_info  **trees);

void compute_trees_matches(char   *filename_halo_root_in,
                           char   *filename_root_out,
                           int     i_read_stop,
                           int     i_read_start,
                           int     i_read_step,
                           int    *n_files,
                           int   **n_subgroups,
                           int   **n_groups,
                           int     n_search);
void compute_trees_horizontal(char   *filename_halos_root_in,
                              char   *filename_cat_root_in,
                              char   *filename_root_matches,
                              char   *filename_root_out,
                              double *a_list,
                              int     i_read_start,
                              int     i_read_stop,
                              int     i_read_step,
                              int     n_search,
                              int    *flag_clean);
void compute_trees_vertical(char *filename_root_out,
                            int   i_read_start,
                            int   i_read_stop,
                            int   i_read_step,
                            int   n_files_groups,
                            int   n_files_subgroups,
                            int  *flag_clean);
void finalize_trees_vertical(tree_info **trees,
                             int        *n_halos_tree,
                             int         n_trees,
                             int         n_snaps,
                             int         progenitor_mode);
void write_trees_vertical(tree_info **trees,
                          int        *n_halos_tree_local,
                          int         n_trees_local,
                          int        *tree_lo_file,
                          int        *tree_hi_file,
                          int        *n_halos_file,
                          int         n_files,
			  char       *filename_root_out,
                          char       *group_text_prefix);
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
#endif

