#ifndef GBPTREES_BUILD_AWAKE
#define GBPTREES_BUILD_AWAKE
#include <gbpHalos.h>
#include <gbpCosmo.h>

// This defines the minimum effective fraction of
//    mass needed to be considered a good match
#define F_GOODNESS_OF_MATCH   0.0
//#define F_GOODNESS_OF_MATCH   0.05

// This defines the minimum and minimum fraction of the 
//    max score needed to be considered a good match
#define F_MAX_MATCH_SCORE_MIN 0.0
#define MIN_MATCH_SCORE       1.0

#define K_MATCH_SUBGROUPS 0
#define K_MATCH_GROUPS    1

// Tree finalization and reading modes
#define TREE_READ_DEFAULT                      0
#define TREE_SUBSTRUCTURE_ORDER_DEFAULT        TTTP00
#define TREE_PROGENITOR_ORDER_DELUCIA          TTTP01
#define TREE_PROGENITOR_ORDER_N_PARTICLES      TTTP02
#define TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK TTTP03
#define TREE_READ_EXTENDED_POINTERS            TTTP04
#define TREE_READ_HEADER_ONLY                  TTTP05
#define TREE_MODE_REFERENCE                    TTTP06
#define TREE_PROGENITOR_ORDER_DEFAULT          TREE_PROGENITOR_ORDER_N_PARTICLES_PEAK
#define TREE_MODE_DEFAULT                      (TREE_SUBSTRUCTURE_ORDER_DEFAULT|TREE_PROGENITOR_ORDER_DEFAULT)

// If any of these are changed, don't forget to modify parse_match_type.c (TTTPXX means "two-to-the-power-XX")
#define TREE_CASE_NO_PROGENITORS               TTTP00  // Set for halos that have no progenitors.
#define TREE_CASE_MAIN_PROGENITOR              TTTP01  // Set for the progenitor with the highest match score. 
#define TREE_CASE_MOST_MASSIVE                 TTTP02  // Marks the most massive substructure.
#define TREE_CASE_DOMINANT                     TTTP03  // Marks the dominant     substructure.
#define TREE_CASE_REMNANT                      TTTP04  // Set for halos with more than one progenitor.
#define TREE_CASE_MERGER_PRIMARY               TTTP05  // Set when a halo is deemed to be the primary   progenitor of a merger
#define TREE_CASE_MERGER                       TTTP06  // Set when a halo is deemed to be the secondary progenitor of a merger
#define TREE_CASE_STRAYED                      TTTP07  // Set for halos for which a descendant was not found
#define TREE_CASE_DROPPED                      TTTP08  // Set if file_offset>1 and TREE_CASE_MATCHED_TO_BRIDGE is not set
#define TREE_CASE_BRIDGED                      TTTP09  // Set for halos with multiple unique back-matches from halos with unique IDs
#define TREE_CASE_EMERGED                      TTTP10  // Set when a match is made identifying this halo as emerged
#define TREE_CASE_FRAGMENTED_NEW               TTTP11  // Set for halos that have been marked TREE_CASE_EMERGED_CANDIDATE but not TREE_CASE_EMERGED
                                                       //    (unless it's the backmatch with the most massive descendant; that halo is considered
                                                       //     to be the source of any fragmented halos)
#define TREE_CASE_FRAGMENTED_STRAYED           TTTP12  // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                       //    decendant_id!=a valid id (ie they are not a progenitor of anything). 
#define TREE_CASE_FRAGMENTED_RETURNED          TTTP13  // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                       //    decendant_id==the id of the halo they are emerged from. 
#define TREE_CASE_FRAGMENTED_EXCHANGED         TTTP14  // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                       //    decendant_id!=the id of the halo they are emerged but is nevertheless valid 
                                                       //    (ie. they are still a progenitor of something). 
#define TREE_CASE_EMERGED_CANDIDATE            TTTP15  // Set when a halo is identified as a unique back-match to a halo marked TREE_CASE_BRIDGED 
                                                       //    and is not identified as the BRIDGE's main descendant
#define TREE_CASE_MATCHED_TO_EMERGED           TTTP16  // Set when a halo is matched to an emerged halo
#define TREE_CASE_2WAY_MATCH                   TTTP17  // Set when the match between a halo and it's descendant is mutual
#define TREE_CASE_GHOST                        TTTP18  // Marks ghost halos in ghost-populated trees
#define TREE_CASE_GHOST_NULL                   TTTP19  // Marks a ghost halo where a subgroup is it's own group.
                                                       //    This is a default behaviour that occurs when a group is strayed but one of 
                                                       //    it's subgroups isn't.
#define TREE_CASE_REINIT_DOMINANT              TTTP20  // For internal use.  This should never be seen in the output.
#define TREE_CASE_UNPROCESSED                  TTTP21  // For internal use.  This should never be seen in the output.
#define TREE_CASE_INVALID                      TTTP22  // For internal use.  This should never be seen in the output.

#ifdef _MAIN
   int   n_tree_case_flag_list=22;
   int   tree_case_flag_list[]={
                  TREE_CASE_2WAY_MATCH,
                  TREE_CASE_NO_PROGENITORS,
                  TREE_CASE_FRAGMENTED_NEW,
                  TREE_CASE_FRAGMENTED_STRAYED,
                  TREE_CASE_FRAGMENTED_RETURNED,
                  TREE_CASE_FRAGMENTED_EXCHANGED,
                  TREE_CASE_MAIN_PROGENITOR,
                  TREE_CASE_MOST_MASSIVE,
                  TREE_CASE_DOMINANT,
                  TREE_CASE_MERGER_PRIMARY,
                  TREE_CASE_MERGER,
                  TREE_CASE_BRIDGED,
                  TREE_CASE_MATCHED_TO_EMERGED,
                  TREE_CASE_REMNANT,
                  TREE_CASE_DROPPED,
                  TREE_CASE_EMERGED,
                  TREE_CASE_EMERGED_CANDIDATE,
                  TREE_CASE_GHOST,
                  TREE_CASE_GHOST_NULL,
                  TREE_CASE_STRAYED,
                  TREE_CASE_UNPROCESSED,
                  TREE_CASE_INVALID};
   const char *tree_case_flag_list_text[]={
                        "2WAY",
                        "NO_PROGENITORS",
                        "FRAGMENTED_NEW",
                        "FRAGMENTED_STRAYED",
                        "FRAGMENTED_RETURNED",
                        "FRAGMENTED_EXCHANGED",
                        "MAIN_PROGENITOR",
                        "MOST_MASSIVE",
                        "DOMINANT",
                        "MERGER_PRIMARY",
                        "MERGER",
                        "BRIDGED",
                        "MATCHED_TO_EMERGED",
                        "REMNANT",
                        "DROPPED",
                        "EMERGED",
                        "EMERGED_CANDIDATE",
                        "GHOST",
                        "GHOST_NULL",
                        "STRAYED",
                        "UNPROCESSED",
                        "INVALID"};
#else
   extern int n_tree_case_flag_list;
   extern int tree_case_flag_list[];
   extern const char *tree_case_flag_list_text[];
#endif

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
#define READ_TREES_CATALOGS_SAGE        4
#define READ_TREES_CATALOGS_BOTH        (READ_TREES_CATALOGS_GROUPS|READ_TREES_CATALOGS_SUBGROUPS)
#define READ_TREES_CATALOGS_ALL         (READ_TREES_CATALOGS_BOTH)|READ_TREES_CATALOGS_PROFILES
#define READ_TREES_CATALOGS_DEFAULT     READ_TREES_CATALOGS_BOTH

#define READ_TREES_DATA_SUBGROUPS   0
#define READ_TREES_DATA_GROUPS      1
#define READ_TREES_DATA_BOTH        (READ_TREES_DATA_GROUPS|READ_TREES_DATA_SUBGROUPS)

#define READ_TREES_MATCH_SCORES_SUBGROUPS   0
#define READ_TREES_MATCH_SCORES_GROUPS      1
#define READ_TREES_MATCH_SCORES_ALL         READ_TREES_MATCH_SCORES_GROUPS|READ_TREES_MATCH_SCORES_SUBGROUPS
#define READ_TREES_MATCH_SCORES_DEFAULT     READ_TREES_MATCH_SCORES_SUBGROUPS

#define READ_TREES_POINTERS_BRIDGE_FOREMATCH 1
#define READ_TREES_POINTERS_BRIDGE_BACKMATCH 2

// Data structures for horizontal tree construction
typedef struct tree_horizontal_stats_info tree_horizontal_stats_info;
struct tree_horizontal_stats_info {
   int n_halos;
   int n_mergers;
   int n_strayed;
   int n_dropped;
   int n_bridged;
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
   int max_emerged_size;
   int max_fragmented_strayed_size;
   int max_fragmented_returned_size;
   int max_fragmented_exchanged_size;
   int max_emerged_progenitor_size;
   int max_id;
};

typedef struct tree_horizontal_info tree_horizontal_info;
typedef struct match_info match_info;
struct match_info{
  tree_horizontal_info *halo;
  float                 score;
  int                   flag_two_way;
};
typedef struct back_match_info back_match_info;
struct back_match_info{
  tree_horizontal_info *halo;
  float                 score;
  int                   file;
};

struct tree_horizontal_info{
  int              id;                             // This halo's id
  int              main_progenitor_id;             // This halo's main progenitor id
  int              tree_id;                        // This halo's tree id
  int              type;                           // A bit-wise switch characterising this halo's matching
  int              file;                           // This halo's snapshot index (ie. 0->n_snaps_used_in_trees-1)
  int              snap;                           // This halo's snapshot number
  int              n_particles;                    // Number of particles in this halo
  int              n_particles_parent;             // Number of particles in this halo's parent halo
  int              n_particles_largest_descendant; // Number of particles in this halo's largest-ever descendant
  int              n_back_matches;                 // The number of halos back-matched to this halo (may contain it's descendant as well)
  int              n_progenitors;                  // The number of progenitors pointing to this halo
  int              index;                          // This halo's index in the halo catalog
  match_info       first_progenitor;               // Pointer to this halo's first progenitor
  match_info       last_progenitor;                // Pointer to this halo's last  progenitor
  match_info       next_progenitor;                // Pointer to this halo's next  progenitor
  back_match_info *back_matches;                   // Contains the pointer information for all of the back-matches to this halo
  match_info       forematch_first;                // Pointer to the first bridged halo matched to.
  match_info       forematch_default;              // Pointer to the default halo matched to.  Starts with an initial match but may change
                                                   //    if we manage to match to one-or-more emerged halo(s).  When we finish scanning, 
                                                   //    this becomes the match.
  match_info       forematch_best;                 // The best match to this halo generated in 'construct_progenitors'.  Used
                                                   //    to ensure that a bridged halo gets at least one progenitor, rather
                                                   //    than have all matches to it end-up as later matches to emerged halos.
  match_info       bridge_backmatch;               // Pointer to a possible back-matched bridged halo
  match_info       descendant;                     // Contains all the needed pointers to the descendant
};

typedef struct tree_horizontal_extended_info tree_horizontal_extended_info;
struct tree_horizontal_extended_info{
  int id;                        // This halo's id
  int tree_id;                   // This halo's tree id
  int descendant_id;             // This halo's main progenitor id
  int parent_id;                 // This halo's parent's main progenitor id
  int type;                      // A bit-wise switch characterising this halo's matching
  int n_particles;               // Number of particles in this halo
  int n_particles_peak;          // Peak number of particles in this halo's progenitor line
  int n_particles_parent;        // Number of particles in this halo's parent
  int n_particles_desc;          // Number of particles in this halo's descendant
  int n_particles_proj;          // Number of particles in this halo's progenitor
  int score_desc;                // Matching score of this halo to it's descendant
  int score_prog;                // Matching score of this halo to it's progenitor
  int snap_bridge;               // Snapshot of any halo that this halo may be back-matched to
  int file_bridge;               // File index of any halo that this halo may be back-matched to
  int index_bridge;              // Index of any bridge this halo may be back-matched to
  int id_bridge;                 // ID of any bridge this halo may be back-matched to
  int substructure_index;        // This halo's substructure index 
  int descendant_file_offset;    // Pointer to this halo's descendant (ie. 0->n_snaps_used_in_trees-1)
  int descendant_index;          // Pointer to this halo's descendant
  int first_progenitor_file;     // Pointer to this halo's first progenitor
  int first_progenitor_index;    // Pointer to this halo's first progenitor
  int next_progenitor_file;      // Pointer to this halo's next  progenitor
  int next_progenitor_index;     // Pointer to this halo's next  progenitor
};

// Data structures for vertical tree construction
typedef struct tree_vertical_node_info tree_vertical_node_info;
struct tree_vertical_node_info{
  halo_properties_SAGE_info  halo;
  int                        depth_first_index;
  int                        group_id;
  int                        halo_id;
  int                        descendant_id;
  int                        n_progenitors;
  tree_vertical_node_info   *descendant;
  tree_vertical_node_info   *progenitor_first;
  tree_vertical_node_info   *progenitor_next;
  tree_vertical_node_info   *progenitor_last;
  tree_vertical_node_info   *group_halo_first;
  tree_vertical_node_info   *group_halo_next;
  tree_vertical_node_info   *neighbour_halo_next;
  tree_vertical_node_info   *next; // Points to the next halo, in the order they are added to the tree
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
  int             depth_first_index;
  int             n_progenitors;
  int             n_substructures;
  int             snap_tree;
  int             file_index;
  int             neighbour_index;
  int             n_particles;
  int             n_particles_peak;
  float           match_score;
  int             halo_ID;
  int             tree_case;
  // Pointers for the substructure heirarchy
  tree_node_info *parent;
  tree_node_info *substructure_first;   // for substructure in this halo's parent
  tree_node_info *substructure_last;    // for substructure in this halo's parent
  tree_node_info *substructure_next;    // for substructure in this halo's parent
  // Merger tree pointers
  tree_node_info *descendant;
  tree_node_info *progenitor_primary;   // Marks the primary halo for mergers
  tree_node_info *progenitor_first;
  tree_node_info *progenitor_last;
  tree_node_info *progenitor_next;
  // Bulk processing pointers
  tree_node_info *next_neighbour;  // This halo's snapshot
  tree_node_info *next_in_forest;  // This halo's forest
};

typedef struct tree_markers_info tree_markers_info;
struct tree_markers_info{
  tree_node_info *branch_leaf;
  tree_node_info *branch_root;
  tree_node_info *descendant;
  tree_node_info *main_progenitor;
  tree_node_info *first_became_satellite;
  tree_node_info *joined_current_parent;
  tree_node_info *peak_mass;
  tree_node_info *half_peak_mass;
  tree_node_info *merger_33pc_remnant; // Placed at the time of the remnant
  tree_node_info *merger_33pc_host;    // Placed at the point when the secondary's peak mass is reached
  tree_node_info *merger_33pc_merger;  // Placed at the point when the secondary's peak mass is reached
  tree_node_info *merger_10pc_remnant; // Placed at the time of the remnant
  tree_node_info *merger_10pc_host;    // Placed at the point when the secondary's peak mass is reached
  tree_node_info *merger_10pc_merger;  // Placed at the point when the secondary's peak mass is reached
  double          M_peak;
};

typedef struct tree_info tree_info;
struct tree_info{
  int              mode;
  // Filename info
  char             filename_root[MAX_FILENAME_LENGTH];
  char             filename_root_horizontal[MAX_FILENAME_LENGTH];
  char             filename_root_horizontal_trees[MAX_FILENAME_LENGTH];
  char             filename_root_analysis[MAX_FILENAME_LENGTH];
  char             name[MAX_FILENAME_LENGTH];
  // Snapshot info
  int              i_read_start;
  int              i_read_stop;
  int              i_read_step;
  int              i_read_last;
  int              n_snaps;
  int              n_wrap;
  int              n_wrap_lookup;
  int              n_search;
  int             *snap_list;
  double          *a_list;
  double          *z_list;
  double          *t_list;
  // Halo and tree counts etc
  int              n_groups_trees;          // No. in the trees
  int              n_groups_trees_local;    // No. in the trees
  int              n_subgroups_trees;       // No. in the trees
  int              n_subgroups_trees_local; // No. in the trees
  int              n_groups_raw;            // No. in the input catalog files
  int              n_groups_raw_local;      // No. in the input catalog files
  int              n_subgroups_raw;         // No. in the input catalog files
  int              n_subgroups_raw_local;   // No. in the input catalog files
  int              n_groups_max;
  int              n_subgroups_max;
  int              n_progenitors_max;
  int              n_trees_subgroup;
  int              n_trees_group;
  int              max_n_groups_snap_local;
  int              max_n_subgroups_snap_local;
  int              n_groups_snap_alloc_local;
  int              n_subgroups_snap_alloc_local;
  int             *n_groups_snap_local;
  int             *n_subgroups_snap_local;
  int             *n_groups_catalog;
  int             *n_subgroups_catalog;
  // Forest info
  int              n_forests;
  int              n_forests_local;
  int              n_trees_forest_groups_max;
  int              n_trees_forest_subgroups_max;
  int              max_n_groups_forest_local;
  int              max_n_subgroups_forest_local;
  int              forest_lo_group_local;
  int              forest_hi_group_local;
  int              forest_lo_subgroup_local;
  int              forest_hi_subgroup_local;
  int             *n_groups_forest_local;
  int             *n_subgroups_forest_local;
  // Flags
  int              flag_fix_bridges;
  int              flag_compute_fragmented;
  int              flag_compute_ghosts;
  // Pointers
  tree_node_info **first_neighbour_groups;
  tree_node_info **first_neighbour_subgroups;
  tree_node_info **last_neighbour_groups;
  tree_node_info **last_neighbour_subgroups;
  tree_node_info **first_in_forest_groups;
  tree_node_info **first_in_forest_subgroups;
  tree_node_info **last_in_forest_groups;
  tree_node_info **last_in_forest_subgroups;
  // Look-up table stuff for tieing file indices to stored halo structures
  //   (mostly needed just for reading)
  int             **group_indices;
  tree_node_info ***group_array;
  int             **subgroup_indices;
  tree_node_info ***subgroup_array;
  int              *tree2forest_mapping_group;
  int              *tree2forest_mapping_subgroup;
  // An ADaPS structure for holding ancillary data
  ADaPS            *data;
  // Cosmology
  cosmo_info       *cosmo;
  double            box_size;
  // Short-cuts to some common stuff added to data
  float                      **group_match_scores;
  float                      **subgroup_match_scores;
  tree_markers_info          **group_markers;
  tree_markers_info          **subgroup_markers;
  halo_properties_info       **group_properties;
  halo_properties_info       **subgroup_properties;
  halo_properties_SAGE_info  **group_properties_SAGE;
  halo_properties_SAGE_info  **subgroup_properties_SAGE;
  tree_node_info            ***group_backmatch_pointers;
  tree_node_info            ***subgroup_backmatch_pointers;
  tree_node_info            ***group_forematch_pointers;
  tree_node_info            ***subgroup_forematch_pointers;
  float                      **group_descendant_score;
  float                      **subgroup_descendant_score;
  float                      **group_progenitor_score;
  float                      **subgroup_progenitor_score;
  float                      **group_backmatch_score;
  float                      **subgroup_backmatch_score;
  float                      **group_forematch_score;
  float                      **subgroup_forematch_score;
  // Reference trees pointer permitting communication between difference tree versions
  tree_info *trees_reference;
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
int tree_case_flags_text(int match_type,const char *separator_string,char **return_string);
void read_trees(char       *filename_SSiMPL_root,
                char       *filename_halos_version,
                char       *filename_trees_version,
                int         read_mode,
                tree_info **trees);
void read_trees_data(tree_info    *trees,
                     char         *filename_root,
                     int           mode,
                     SID_Datatype  data_type,
                     const char   *name);
void init_trees_data(tree_info    *trees,
                     void       ***rval,
                     size_t        data_size,
                     int           mode,
                     const char   *name,
                     ...);
void init_trees_lookup(tree_info *trees);
void update_trees_lookup(tree_info *trees,int i_file);
int  find_tree_node(tree_info       *trees,
                    int              node_file,
                    int              node_index,
                    int              group_mode,
                    tree_node_info **descendant);
void free_trees_lookup(tree_info *trees);
void free_trees_data(void **tree_data,void *params);
void read_trees_match_scores(tree_info *trees,
                             char      *filename_SSimPL_dir,
                             int        mode);
void read_trees_pointers(tree_info        *trees,
                         const char       *filename_input_dir_horizontal_trees,
                         int               i_file_ptrs,
                         int               i_read_ptrs,
                         int               mode);
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
                  size_t  *match_index,
                  char    *match_flag_two_way,
                  double   f_goodness_of_match);
int check_for_matching_input_files(const char *filename_root_in,int i_read);

float maximum_match_score(double n_particles);
float minimum_match_score(double n_particles);
float match_score_f_goodness(float match_score,int n_particles_in);
int check_goodness_of_match(int n_particles_i,int n_particles_j,float match_score,double f_goodness_of_match);

int check_if_halo_is_descendant(tree_horizontal_info *possible_progenitor,
                                tree_horizontal_info *possible_descendant,
                                int n_search);
int check_if_halo_is_merger(int type);
int check_if_descendant_is_back_matched(tree_horizontal_info *halo,
                                        tree_horizontal_info *halo_to_check);
int check_validity_of_main_progenitor(tree_horizontal_info *descendant_halo,
                                      match_info           *old_MP,
                                      match_info           *new_MP);
int check_validity_of_emerged_match(tree_horizontal_info *halo_i,
                                    back_match_info      *back_match,
                                    char                  match_flag_two_way,
                                    int                   n_search);
int check_validity_of_tree_case_flag(int flag);
int compute_cross_catalog_matches(char   *filename_root_in_1,
                                  char   *filename_root_in_2,
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
                         const char *filename_cat1,
                         const char *filename_cat2,
                         plist_info *plist1,
                         plist_info *plist2,
                         int         k_match,
                         int         mode);
void compute_trees_horizontal_stats(void *halos_in,int n_halos,int n_halos_max,tree_horizontal_stats_info *stats,int flag_write_cases);

int   set_match_id         (match_info *match);
int   set_match_file       (match_info *match);
int   set_match_snapshot   (match_info *match);
int   set_match_type       (match_info *match);
int   set_match_n_particles(match_info *match);
int   set_match_index      (match_info *match);
float set_match_score      (match_info *match);

int   set_back_match_id         (back_match_info *back_match);
int   set_back_match_file       (back_match_info *back_match);
int   set_back_match_snapshot   (back_match_info *back_match);
int   set_back_match_type       (back_match_info *back_match);
int   set_back_match_n_particles(back_match_info *back_match);
int   set_back_match_index      (back_match_info *back_match);
float set_back_match_score      (back_match_info *back_match);

void set_n_particles_peak(int type,int n_particles_halo,int *n_particles_peak);

tree_horizontal_extended_info *set_extended_descendant(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int i_file,int n_wrap);
tree_horizontal_extended_info *set_extended_first_progenitor(tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int n_wrap);
tree_horizontal_extended_info *set_extended_next_progenitor (tree_horizontal_extended_info **halos,tree_horizontal_extended_info *halo,int n_wrap);
void check_for_fragmented_halos(int k_match,tree_horizontal_info **groups,int n_groups,
                                int i_write,int j_write,int l_write,int n_wrap);
void add_to_trees_horizontal_stats(tree_horizontal_stats_info *stats,int id,int type,int n_particles);
void init_trees_horizontal_roots(tree_horizontal_info **groups,
                                 tree_horizontal_info **subgroups,
                                 int    *match_id,
                                 float  *match_score,
                                 size_t *match_index,
                                 char   *match_flag_two_way,
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
void init_trees_horizontal_snapshot(tree_horizontal_info *halos,
                                    int                   i_read,
                                    int                   i_file,
                                    int                   n_groups,
                                    int                   n_groups_max,
                                    int                   n_subgroups,
                                    int                   n_subgroups_max,
                                    int                   flag_match_subgroups);
void add_substructure_info(tree_horizontal_info *halos,
                           int                  *n_subgroups_group,
                           int                  *n_particles_groups,
                           int                   n_groups,
                           int                   n_subgroups,
                           int                   flag_match_halos);
void read_halo_sizes(tree_horizontal_info **halos,
                     int     n_halos_i,
                     int    *match_id,
                     float  *match_score,
                     size_t *match_index,
                     int    *n_particles,
                     int     i_file,
                     int     i_read,
                     int     i_read_step,
                     int     n_wrap,
                     int     n_halos_max,
                     char   *filename_root_matches,
                     int     flag_match_subgroups);
void identify_bridges(tree_horizontal_info **halos,
                      tree_horizontal_info  *halos_i,
                      int     n_halos_i,
                      int    *match_id,
                      float  *match_score,
                      size_t *match_index,
                      char   *match_flag_two_way,
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
void identify_emerged_halo_candidates(tree_horizontal_info *halos_i,
                                      int                   n_halos_i,
                                      int                   n_search);
void construct_progenitors(tree_horizontal_info **halos,
                           tree_horizontal_info  *halos_i,
                           int   **n_subgroups_group,
                           int     n_halos_i,
                           int    *match_id,
                           float  *match_score,
                           size_t *match_index,
                           char   *match_flag_two_way,
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
void finalize_bridged_halo_list(tree_horizontal_info *halos_i,
                                int                   n_halos_i,
                                int                   i_file,
                                int                   n_search,
                                int                   n_files);
void set_largest_descendants(tree_horizontal_info *halos_i,
                             int                   n_halos_i);
void finalize_trees_horizontal(int                    n_halos_1_matches,
                               int                    n_halos_i,
                               tree_horizontal_info **halos,
                               tree_horizontal_info  *halos_i,
                               int                    i_file,
                               int                    n_search,
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
void read_forests(const char  *filename_root_in,
                  int   *n_forests_group,
                  int   *n_forests_subgroup,
                  int   *n_forests_group_local,
                  int   *n_forests_subgroup_local,
                  int  **i_forest_group,
                  int  **i_forest_subgroup,
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
void read_trees_final_totals(char *filename_output_dir_horizontal_trees,
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
void propagate_progenitor_info(int         *n_groups,
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
                               char        *filename_output_dir,
                               int          flag_compute_fragmented);
void propagate_dominant_substructures(tree_horizontal_extended_info **groups,   int *n_groups,
                                      tree_horizontal_extended_info **subgroups,int *n_subgroups,
                                      int        **n_subgroups_group,
                                      int          i_read, // tree snapshot index
                                      int          j_read, // actual snapshot index
                                      int          l_read,
                                      int          i_read_step,
                                      int          n_wrap);
void propagate_n_particles_peak(tree_horizontal_extended_info **groups,   int *n_groups,
                                tree_horizontal_extended_info **subgroups,int *n_subgroups,
                                int        **n_subgroups_group,
                                int          i_read, // tree snapshot index
                                int          j_read, // actual snapshot index
                                int          l_read,
                                int          i_read_step,
                                int          n_wrap);
void propagate_bridge_info(tree_horizontal_extended_info **groups,   int *n_groups,
                           tree_horizontal_extended_info **subgroups,int *n_subgroups,
                           int        **n_subgroups_group,
                           int          i_read, // tree snapshot index
                           int          j_read,
                           int          l_read,
                           int          i_read_step,
                           int          n_wrap);
void propagate_fragmented_info(tree_horizontal_extended_info **groups,   int *n_groups,
                               tree_horizontal_extended_info **subgroups,int *n_subgroups,
                               int        **n_subgroups_group,
                               int          i_read, // tree snapshot index
                               int          j_read,
                               int          l_read,
                               int          i_read_step,
                               int          n_wrap);
void propagate_merger_info(tree_horizontal_extended_info **groups,   int *n_groups,
                           tree_horizontal_extended_info **subgroups,int *n_subgroups,
                           int        **n_subgroups_group,
                           int          i_read, // tree snapshot index
                           int          j_read,
                           int          l_read,
                           int          i_read_step,
                           int          n_wrap);
void add_progenitor_to_halo(tree_horizontal_info **halos,
                            int                    i_file,
                            int                    i_halo,
                            int                    j_file,
                            int                    j_halo,
                            float                  score,
                            int                   *max_id,
                            int                    n_wrap);
void compute_trees_horizontal(char   *filename_halos_root_in,
                              char   *filename_cat_root_in,
                              char   *filename_snap_list_in,
                              char   *filename_root_matches,
                              char   *filename_root_out,
                              cosmo_info **cosmo,
                              int     i_read_start,
                              int     i_read_stop,
                              int     i_read_step,
                              int     n_search,
                              int     flag_fix_bridges,
                              int    *flag_clean);
void compute_forests(char *filename_root_out,int n_search_forests);
void compute_trees_vertical(char   *filename_SSimPL_dir,
                            char   *filename_halo_version_root,
                            char   *filename_trees_name,
                            double  box_size,
                            int     n_dim_files);
void finalize_trees(tree_info *trees,
                    int        progenitor_mode);
void finalize_trees_vertical(tree_info *trees);
void read_trees_run_parameters(const char *filename_root_out,
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
void write_trees_vertical(tree_info     *trees,
                          double         box_size,
                          int            grid_size,
                          const char    *filename_root_out);
void write_a_list(const char *filename_snap_list_in,
                  const char *filename_root_out,
                  int         i_read_start,
                  int         i_read_stop,
                  int         i_read_step);
void read_a_list(const char  *filename_root_in,
                 double     **a_list_in,
                 int         *n_a_list);
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
void init_trees_read(const char  *filename_SSimPL_dir,
                     const char  *filename_trees_version,
                     int          mode,
                     tree_info  **tree);
void free_trees(tree_info **tree);
int add_node_to_trees(tree_info        *trees,
                      int               i_forest,
                      int               tree_case,
                      int               n_particles,
                      int               halo_ID,
                      int               halo_snap,
                      int               halo_index,
                      int               descendant_snap,
                      int               descendant_index,
                      tree_node_info   *group_node,
                      tree_node_info  **new_node);
void init_trees_vertical(int n_snaps,tree_vertical_info **tree);
void free_trees_vertical(tree_vertical_info **tree);
int  add_node_to_vertical_tree(tree_vertical_info        *tree,
                               int                        match_type,
                               int                        halo_id,
                               int                        group_id,
                               int                        descendant_id,
                               int                        halo_snap,
                               int                        descendant_snap,
                               halo_properties_SAGE_info *properties);

void assign_group_order(tree_info *tree,int mode);
void assign_progenitor_order(tree_info *tree,int mode);
void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode);
void assign_depth_first_index_vertical_recursive(tree_node_info *tree,int *depth_first_index);
void assign_unique_vertical_tree_ids(tree_info *trees,tree_node_info *tree_node);
void compute_progenitor_score_recursive(tree_node_info *tree,int *M_i,int mode);
void compute_substructure_order_recursive(tree_node_info *parent,int *score_parent,int mode);
void compute_progenitor_order_recursive(tree_info *trees,tree_node_info *descendant,int *score_descendant,int mode);
int  construct_unique_vertical_tree_id(tree_node_info *tree_node);

int  construct_unique_tree_id(tree_node_info *tree_node);
void assign_depth_first_index_recursive(tree_node_info *tree,int *depth_first_index);
void assign_unique_tree_ids_recursive(tree_node_info *tree_node,int i_tree);
void assign_group_subgroup_order_vertical(tree_vertical_info *tree,int i_snap,int mode);
void assign_progenitor_order_vertical_recursive(tree_vertical_node_info *tree,int *M_i,int mode);
void compute_progenitor_score_vertical_recursive(tree_vertical_node_info *tree,int *M_i,int mode);

#ifdef __cplusplus
}
#endif

#endif

