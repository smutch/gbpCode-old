
#define TREE_PROGENITOR_ORDER_DEFAULT 0
#define TREE_PROGENITOR_ORDER_DELUCIA 2

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

void compute_trees_horizontal(char *filename_halos_root_in,
                              char *filename_cat_root_in,
                              char *filename_root_out,
                              int   i_read_start,
                              int   i_read_stop,
                              int   i_read_step,
                              int   n_search,
                              int  *flag_clean);
void compute_trees_vertical(char *filename_root_out,
                            int   i_read_start,
                            int   i_read_stop,
                            int   i_read_step,
                            int   n_search,
                            int   n_files_groups,
                            int   n_files_subgroups,
                            int  *flag_clean);
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

