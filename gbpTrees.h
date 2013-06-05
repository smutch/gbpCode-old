#ifndef GBPTREES_AWAKE
#define GBPTREES_AWAKE
#include <gbpHalos.h>
#include <gbpCosmo.h>

#define TREE_PROGENITOR_ORDER_DEFAULT 0
#define TREE_PROGENITOR_ORDER_DELUCIA 2

// If any of these are changed, don't forget to modify parse_match_type.c
#define TREE_CASE_SIMPLE                        1       // Set when a halo has a file_offset=1 
#define TREE_CASE_MAIN_PROGENITOR               2       // Set for the progenitor with the highest match score. (propagated for ghosts)
#define TREE_CASE_MERGER                        4       // Set when new IDs are created (ie. last point the halo was seen).
                                                        //    Set only for the last ghost in ghost-populated trees for mergers w/ offset>1.
#define TREE_CASE_DROPPED                       8       // Set if file_offset>1 and TREE_CASE_MATCHED_TO_BRIDGE is not set
#define TREE_CASE_STRAYED                       16      // Set for halos for which a descendant was not found
#define TREE_CASE_SPUTTERED                     32      // Set for halos whose descendant was not given a valid ID. (propagated for ghosts)
#define TREE_CASE_BRIDGED                       64      // Set for halos with multiple back-matches from halos with unique IDs
#define TREE_CASE_EMERGED_CANDIDATE             128     // Set when a halo is identified as a unique back-match to a halo marked TREE_CASE_BRIDGED 
                                                        //    and is not identified as the BRIDGE's main descendant
#define TREE_CASE_FOUND                         256     // Set when a halo has a progenitor with a file_offset>1
#define TREE_CASE_NO_PROGENITORS                512     // Set for halos that have no progenitors.
#define TREE_CASE_FRAGMENTED_LOST               1024    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                        //    decendant_id!=a valid id (ie they are not a progenitor of anything). (propagated for ghosts)
#define TREE_CASE_FRAGMENTED_RETURNED           2048    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                        //    decendant_id==the id of the halo they are emerged from. (propagated for ghosts)
#define TREE_CASE_FRAGMENTED_EXCHANGED          4096    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                        //    decendant_id!=the id of the halo they are emerged but is nevertheless valid 
                                                        //    (ie. they are still a progenitor of something). (propagated for ghosts)
#define TREE_CASE_MATCHED_TO_BRIDGE             8192    // Set when a halo is matched to one with TREE_CASE_BRIDGED set
#define TREE_CASE_BRIDGE_DEFAULT                16384   // Set when a halo matched to a bridge is not matched to any emerged candidate halos
#define TREE_CASE_GHOST                         32768   // Marks ghost halos in ghost-populated trees
#define TREE_CASE_GHOST_NAKED                   65536   // Marks a ghost halo where a subgroup is it's own group (product of a default action).
#define TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED 131072  // For internal use.  This should never be seen in the output.
#define TREE_CASE_BRIDGE_FINALIZE               262144  // For internal use.  This should never be seen in the output.
#define TREE_CASE_UNPROCESSED                   524288  // For internal use.  This should never be seen in the output.
#define TREE_CASE_INVALID                       1048576 // For internal use.  This should never be seen in the output.
#define TREE_CASE_EMERGED                       (TREE_CASE_EMERGED_CANDIDATE+TREE_CASE_FOUND)
#define TREE_CASE_FRAGMENTED_NEW                (TREE_CASE_EMERGED_CANDIDATE+TREE_CASE_NO_PROGENITORS)

#define TREE_HORIZONTAL_READ_DEFAULT   0
#define TREE_HORIZONTAL_READ_EXTENDED  1
#define TREE_HORIZONTAL_STORE_EXTENDED 2
#define TREE_HORIZONTAL_STORE_GHOSTS   4

#define TREE_HORIZONTAL_WRITE_DEFAULT   0
#define TREE_HORIZONTAL_WRITE_ALLCASES  1
#define TREE_HORIZONTAL_WRITE_NOCASES   2
#define TREE_HORIZONTAL_WRITE_EXTENDED  4
#define TREE_HORIZONTAL_WRITE_GHOSTS    8

// Data structures for horizontal tree construction
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
   int n_fragmented_lost;
   int n_fragmented_returned;
   int n_fragmented_exchanged;
   int n_emerged_progenitors;
   int n_invalid;
   int n_unprocessed;
   int max_strayed_size;
   int max_sputtered_size;
   int max_dropped_size;
   int max_bridged_size;
   int max_bridge_progenitor_size;
   int max_emerged_size;
   int max_fragmented_lost_size;
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

typedef struct tree_horizontal_read_info tree_horizontal_read_info;
struct tree_horizontal_read_info{
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
  int          file_bridge;        // File index of any halo that this halo may be matched to
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
typedef struct tree_node_info tree_node_info;
struct tree_node_info{
  halo_info       halo;
  int             depth_first_index;
  int             group_id;
  int             halo_id;
  int             descendant_id;
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

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
void read_trees(char             *filename_trees_root,
                int               i_file_start,
                int               i_file_stop,
                int               mode,
                tree_node_info  **trees);
void read_AHF_for_trees(char       *filename_root,
                        int         i_file,
                        plist_info *plist,
                        char       *catalog_name,
                        int         mode);
void read_matches(char    *filename_root_matches,
                  int      i_read,
                  int      j_read,
                  int      mode,
                  int     *n_groups_i,
                  int     *n_groups_j,
                  int     *n_particles_i,
                  int     *n_particles_j,
                  int     *n_sub_group_i,
                  int     *n_sub_group_j,
                  int     *match_ids,
                  float   *match_score,
                  size_t  *match_index);
void compute_trees_matches(char   *filename_halo_root_in,
                           char   *filename_root_out,
                           int     i_read_stop,
                           int     i_read_start,
                           int     i_read_step,
                           int    *n_files,
                           int   **n_subgroups,
                           int   **n_groups,
                           int     n_search);
void write_match_results(char       *filename_out_dir,
                         char       *filename_out_root,
                         int         i_read,
                         int         j_read,
                         plist_info *plist1,
                         plist_info *plist2,
                         int         k_match);
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
                               int i_write,int j_write,int n_wrap);
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
                                    int                   *max_id);
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
                            int     mode);
void write_trees_horizontal_emerged_candidates(int                   i_read,
                                               int                   n_halos_i,
                                               tree_horizontal_info *halos_i,
                                               char                 *group_text_prefix,
                                               char                 *filename_output_dir,
                                               int                   flag_start_new_file);
void count_ghosts(int  *n_groups_in,    int *n_group_ghosts,
                  int  *n_subgroups_in, int *n_subgroup_ghosts,
                  char *filename_output_dir,
                  int   i_file,int i_read);
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
void propagate_fragmented_halos(tree_horizontal_read_info **groups,   int *n_groups,
                                tree_horizontal_read_info **subgroups,int *n_subgroups,
                                int        **n_subgroups_group,
                                int          i_read, // tree snapshot index
                                int          j_read,
                                int          l_read,
                                int          n_wrap);
void set_halo_and_descendant(tree_horizontal_info **halos,
                             int                    i_file,
                             int                    i_halo,
                             int                    j_file,
                             int                    j_halo,
                             float                  score,
                             int                   *max_id,
                             int                    n_wrap);
void create_ghosts(tree_horizontal_ghost_group_info    **groups,
                   tree_horizontal_ghost_subgroup_info **subgroups,
                   int        *n_groups,
                   int        *n_subgroups,
                   int       **n_subgroups_group,
                   int        *n_group_ghosts_used,
                   int        *n_subgroup_ghosts_used,
                   int         i_read,
                   int         l_read,
                   int         j_file,
                   int         i_file_start,
                   int         n_search,
                   int         n_wrap,
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
void compute_trees_vertical(char *filename_root_out,
                            char *filename_cat_root_in,
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
			  const char *filename_root_out,
                          const char *group_text_prefix);
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

void init_tree(int n_snaps,tree_info **tree);
void free_tree(tree_info **tree);
void add_node_to_tree(tree_info  *tree,
                      int         match_type,
                      int         halo_id,
                      int         group_id,
                      int         descendant_id,
                      int         halo_snap,
                      int         descendant_snap,
                      halo_info  *properties);
int  construct_unique_id(tree_node_info *tree_node,int tree_number);
void compute_halo_score_recursive(tree_node_info *tree,int *M_i,int mode);
void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode);
void assign_group_halo_order(tree_info *tree,int i_snap,int mode);
void assign_depth_first_index_recursive(tree_node_info *tree,int *depth_first_index);
void assign_unique_ids_recursive(tree_node_info *tree_node,int i_tree);

#ifdef __cplusplus
}
#endif

#endif

