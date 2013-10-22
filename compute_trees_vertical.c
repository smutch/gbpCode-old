#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

#define K_MATCH_SUBGROUPS 0
#define K_MATCH_GROUPS    1

void init_tree(int n_snaps,tree_info **tree){
  int i_search;
  (*tree)                     =(tree_info       *)SID_malloc(sizeof(tree_info));
  (*tree)->n_neighbours       =(int             *)SID_malloc(sizeof(int)*n_snaps);
  (*tree)->neighbour_halos    =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_snaps);
  (*tree)->neighbour_halo_last=(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_snaps);
  for(i_search=0;i_search<n_snaps;i_search++){
    (*tree)->n_neighbours[i_search]       =0;
    (*tree)->neighbour_halos[i_search]    =NULL;
    (*tree)->neighbour_halo_last[i_search]=NULL;
  }
  (*tree)->root     =NULL;
  (*tree)->last_leaf=NULL;
}

void free_tree(tree_info **tree){
  tree_node_info *next_node;
  tree_node_info *last_node;
  // Free nodes
  next_node=(*tree)->root;
  while(next_node!=NULL){
    last_node=next_node;
    next_node=last_node->next;
    SID_free(SID_FARG last_node);
  }
  // Free neighbour list arrays
  SID_free(SID_FARG (*tree)->n_neighbours);
  SID_free(SID_FARG (*tree)->neighbour_halos);
  SID_free(SID_FARG (*tree)->neighbour_halo_last);
  SID_free(SID_FARG (*tree));
}

int add_node_to_tree(tree_info  *tree,
                     int         match_type,
                     int         halo_id,
                     int         group_id,
                     int         descendant_id,
                     int         halo_snap,
                     int         descendant_snap,                      
                     halo_info  *properties){
  tree_node_info *new_node;
  tree_node_info *last_node;
  tree_node_info *next_node;
  tree_node_info *last_progenitor;
  tree_node_info *next_progenitor;
  tree_node_info *descendant_halo_list;
  tree_node_info *group_halo_list;

  // Create new node
  new_node=(tree_node_info *)SID_malloc(sizeof(tree_node_info));

  // Copy halo properties into new node
  memcpy(&(new_node->halo),properties,sizeof(halo_info));
  
  // Set ids and pointer defaults for the new node
  new_node->depth_first_index    =            -1; // Default; back-filled later
  new_node->n_progenitors        =             0;
  new_node->group_id             =      group_id;
  new_node->halo_id              =       halo_id;
  new_node->descendant_id        = descendant_id;
  new_node->descendant           =          NULL; // Set below
  new_node->progenitor_first     =          NULL; // Set below
  new_node->progenitor_next      =          NULL; // Set below
  new_node->progenitor_last      =          NULL; // Set below
  new_node->group_halo_first     =      new_node; // Set below if group_id>=0
  new_node->group_halo_next      =          NULL; // Set below
  new_node->neighbour_halo_next  =          NULL; // Set below
  new_node->next                 =          NULL;

  // Find the halo's descendant and set various pointers
  //   All used halos must be added to the neighbour list (done below) for this to work
  int flag_found_group=FALSE;
  if(halo_snap!=descendant_snap && descendant_snap>=0){
    descendant_halo_list=tree->neighbour_halos[descendant_snap];
    if(descendant_halo_list!=NULL && descendant_id>=0){
      next_node=descendant_halo_list;
      while(next_node!=NULL){
        if(next_node->halo_id==descendant_id){
          // Set descendant pointer and it's last progenitor
          new_node->descendant=next_node;
          new_node->descendant->n_progenitors++;
          // If this is a descendant's first progenitor, then initialize list ...
          if(new_node->descendant->progenitor_first==NULL){
            new_node->descendant->progenitor_first=new_node;
            new_node->descendant->progenitor_last =new_node;
          }
          // ... else add the new halo to the end of the list
          else{
            next_progenitor=new_node->descendant->progenitor_first;
            while(next_progenitor!=NULL){
              next_progenitor->progenitor_last=new_node;
              last_progenitor                 =next_progenitor;
              next_progenitor                 =next_progenitor->progenitor_next;
            }
            last_progenitor->progenitor_next=new_node;
          }
          break;
        }
        else
          next_node=next_node->neighbour_halo_next;
      }
      if(new_node->descendant==NULL)
        SID_trap_error("Could not find the descendant (%d;snap=%d) of a halo (%d;snap=%d)!",ERROR_LOGIC,halo_id,halo_snap,descendant_id,descendant_snap);
    }
  }
  else if(descendant_snap==halo_snap)
    SID_trap_error("This halo (%d;snap=%d) has the same snap as its descenant (%d;snap=%d)!",ERROR_LOGIC,halo_id,halo_snap,descendant_id,descendant_snap);

  // Find the new halo's group and then set the group pointers
  //   If group_id<0, then this halo is it's own group and we 
  //   can keep defaults and skip this
  group_halo_list=tree->neighbour_halos[halo_snap];
  if(group_halo_list!=NULL && group_id>=0){
    // Scan all the groups in this snapshot
    next_node=group_halo_list;
    while(next_node!=NULL){
      // If we've found a halo belonging to the same group as the new node...
      if(next_node->group_id==group_id){
        // ... then set the new node's group_halo_first pointer ...
        new_node->group_halo_first=next_node->group_halo_first;
        // ... and set the group_halo_next of the last halo that was added to the group ...
        last_node=next_node->group_halo_first;
        next_node=last_node->group_halo_next;
        while(next_node!=NULL){ // we exit the previous while automatically this way
          last_node=next_node;
          next_node=last_node->group_halo_next;
        }
        last_node->group_halo_next=new_node;
        flag_found_group=TRUE;
      }
      // ... else, keep searching
      else
        next_node=next_node->neighbour_halo_next;
    }
  }

  // Set tree neighbour pointers (needed for finding groups and descendants)
  if(tree->neighbour_halos[halo_snap]==NULL)
    tree->neighbour_halos[halo_snap]=new_node;
  else
    tree->neighbour_halo_last[halo_snap]->neighbour_halo_next=new_node;
  tree->neighbour_halo_last[halo_snap]=new_node;
  tree->n_neighbours[halo_snap]++;

  // Set the tree's root and leaf pointers (needed later for whole-tree processing)
  if(tree->root==NULL)
    tree->root=new_node;
  if(tree->last_leaf!=NULL)
    tree->last_leaf->next=new_node;
  tree->last_leaf=new_node;

  return(flag_found_group);
}

int construct_unique_id(tree_node_info *tree_node,int tree_number){
  if(tree_node!=NULL){
    //return(tree_number*1000000+tree_node->depth_first_index);
    return(tree_node->depth_first_index);
  }
  else
    return(-1);    
}

void compute_halo_score_recursive(tree_node_info *tree,int *M_i,int mode){
  tree_node_info  *current;
  int              i_progenitor;
  int              M_iN,N_i,max_M_iN; // Defined in Section 2 of De Lucia and Blaizot (2006)

  N_i     =tree->halo.n_particles;
  max_M_iN=0;
  if(tree->n_progenitors>=1){
    current=tree->progenitor_first;
    i_progenitor=0;
    while(current!=NULL){
      M_iN=0;
      compute_halo_score_recursive(current,&M_iN,mode);
      if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
        max_M_iN=MAX(max_M_iN,M_iN);
      i_progenitor++;
      current=current->progenitor_next;
    }
    if(i_progenitor!=tree->n_progenitors)
      SID_trap_error("There is a progenitor problem in compute_halo_score_recursive! (%d!=%d)",ERROR_LOGIC,i_progenitor!=tree->n_progenitors);
  }

  // Add this progenitor's score to the descendant's sum (see De Lucia and Blaizot (2006))
  (*M_i)=N_i+max_M_iN;
}

void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode){
  tree_node_info  *first_new;
  tree_node_info  *last_new;
  tree_node_info  *current;
  tree_node_info **progenitors;
  int              i_progenitor;
  size_t          *M_iN_index;
  size_t          *M_iN_rank; 
  int             *M_iN,N_i,max_M_iN; // Defined in Section 2 of De Lucia and Blaizot (2006)

  N_i     =tree->halo.n_particles;
  max_M_iN=0;
  if(tree->n_progenitors>=1){
    // Sum all progenitor contributions to score
    M_iN        =(int             *)SID_malloc(sizeof(int)*tree->n_progenitors);
    progenitors =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*tree->n_progenitors);
    i_progenitor=0;
    current     =tree->progenitor_first;
    while(current!=NULL){
      M_iN[i_progenitor]=0;
      assign_progenitor_order_recursive(current,&(M_iN[i_progenitor]),mode);
      progenitors[i_progenitor]=current;
      if(check_mode_for_flag(mode,TREE_PROGENITOR_ORDER_DELUCIA))
        max_M_iN=MAX(max_M_iN,M_iN[i_progenitor]);
      i_progenitor++;
      current=current->progenitor_next;
    }
    if(i_progenitor!=tree->n_progenitors)
      SID_trap_error("There is a progenitor problem in assign_progenitor_order_recursive! (%d!=%d)",ERROR_LOGIC,i_progenitor!=tree->n_progenitors);

    // Assign progenitors by (descending) order of their score
    //   (note: sorting the sort indicies gives each progenitor's ranking in the sort)
    merge_sort(M_iN,      (size_t)tree->n_progenitors,&M_iN_index,SID_INT,   SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    merge_sort(M_iN_index,(size_t)tree->n_progenitors,&M_iN_rank, SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    first_new=progenitors[M_iN_index[tree->n_progenitors-1]]; 
    last_new =progenitors[M_iN_index[0]]; 
    tree->progenitor_first=first_new;
    tree->progenitor_last =last_new;
    for(i_progenitor=0;i_progenitor<tree->n_progenitors;i_progenitor++){
      if(progenitors[i_progenitor]!=last_new)
        progenitors[i_progenitor]->progenitor_next=progenitors[M_iN_index[M_iN_rank[i_progenitor]-1]]; 
      else
        progenitors[i_progenitor]->progenitor_next=NULL;
    }

    SID_free(SID_FARG progenitors);
    SID_free(SID_FARG M_iN);
    SID_free(SID_FARG M_iN_index);
    SID_free(SID_FARG M_iN_rank);
  }

  // Add this progenitor's score to the descendant's sum (see De Lucia and Blaizot (2006))
  (*M_i)=N_i+max_M_iN;
}

void assign_group_halo_order(tree_info *tree,int i_snap,int mode){
  int              n_neighbours;
  int              i_halo,j_halo,k_halo;
  int             *group_ids;
  int              largest_group;
  int              n_in_group;
  tree_node_info **neighbours;
  tree_node_info  *first_neighbour;
  tree_node_info  *current;
  tree_node_info  *new_first;
  tree_node_info  *new_last;
  int             *halo_score;
  size_t          *group_ids_index;
  size_t          *halo_score_index;
  size_t          *halo_score_rank;

  n_neighbours=tree->n_neighbours[i_snap];
  if(n_neighbours>0){
    // Initialize some temporary arrays
    group_ids      =(int             *)SID_malloc(sizeof(int)*n_neighbours);
    neighbours     =(tree_node_info **)SID_malloc(sizeof(tree_node_info)*n_neighbours);
    first_neighbour=tree->neighbour_halos[i_snap];

    // Sort halos by group_id
    i_halo =0;
    current=first_neighbour;
    while(current!=NULL){
      if(i_halo>=n_neighbours)
        SID_trap_error("There's a problem with the number of neighbours in assign_group_halo_order (%d>=%d)!",ERROR_LOGIC,i_halo,n_neighbours);
      neighbours[i_halo]=current;
      group_ids[i_halo] =current->group_id;
      current           =current->neighbour_halo_next;
      i_halo++;
    }
    if(i_halo!=n_neighbours)
      SID_trap_error("There's a problem with the number of neighbours in assign_group_halo_order (%d!=%d)!",ERROR_LOGIC,i_halo,n_neighbours);
    merge_sort(group_ids,(size_t)n_neighbours,&group_ids_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

    // Find highest group score and initialize a temporary array
    largest_group=0;
    i_halo       =0;
    while(i_halo<n_neighbours){
      n_in_group=0;
      while(group_ids[group_ids_index[i_halo+n_in_group]]==group_ids[group_ids_index[i_halo]]){ 
        n_in_group++;
        if((i_halo+n_in_group)>=n_neighbours)
          break;
      }
      largest_group=MAX(largest_group,n_in_group);
      i_halo      +=n_in_group;
    }

    // Scan neighbour list (now ordered by ID) and order by score within each group
    if(largest_group>1){
      halo_score=(int *)SID_malloc(sizeof(int)*largest_group);
      i_halo   =0;
      while(i_halo<n_neighbours){
        // Create a list of halo scores for each group
        n_in_group=0;
        while(group_ids[group_ids_index[i_halo+n_in_group]]==group_ids[group_ids_index[i_halo]]){
          halo_score[n_in_group]=0;
          compute_halo_score_recursive(neighbours[group_ids_index[i_halo+n_in_group]],&(halo_score[n_in_group]),mode);
          n_in_group++;
          if((i_halo+n_in_group)>=n_neighbours)
            break;
        }
        // Correct the ordering within each group ...
        if(n_in_group>1){
          // ... sort halo scores (ascending) ...
          merge_sort(halo_score,      (size_t)n_in_group,&halo_score_index,SID_INT,   SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          merge_sort(halo_score_index,(size_t)n_in_group,&halo_score_rank, SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          // ... set the new pointers ...
          new_first=neighbours[group_ids_index[i_halo+halo_score_index[n_in_group-1]]];
          new_last =neighbours[group_ids_index[i_halo+halo_score_index[0]]];
          for(j_halo=i_halo,k_halo=0;k_halo<n_in_group;j_halo++,k_halo++){
            if(neighbours[group_ids_index[j_halo]]!=new_last)
              neighbours[group_ids_index[j_halo]]->group_halo_next=neighbours[group_ids_index[i_halo+halo_score_index[halo_score_rank[k_halo]-1]]];
            else
              neighbours[group_ids_index[j_halo]]->group_halo_next=NULL;
            neighbours[group_ids_index[j_halo]]->group_halo_first=new_first;
          }
          SID_free(SID_FARG halo_score_index);
          SID_free(SID_FARG halo_score_rank);
        }
        i_halo+=n_in_group;
      }
      SID_free(SID_FARG halo_score);
    }

    SID_free(SID_FARG group_ids);
    SID_free(SID_FARG group_ids_index);
    SID_free(SID_FARG neighbours);
  }
}

void assign_depth_first_index_recursive(tree_node_info *tree,int *depth_first_index){
  tree_node_info *current;

  // Set and increment the tree index
  tree->depth_first_index=(*depth_first_index)++;

  // Walk the tree
  current=tree->progenitor_first;
  while(current!=NULL){
    assign_depth_first_index_recursive(current,depth_first_index);
    current=current->progenitor_next;
  }
}

void assign_unique_ids_recursive(tree_node_info *tree_node,int i_tree){
  tree_node_info *current;
  
  // Set descendant and progenitor ids
  tree_node->halo.descendant      =construct_unique_id(tree_node->descendant,      i_tree);
  tree_node->halo.progenitor_first=construct_unique_id(tree_node->progenitor_first,i_tree);
  tree_node->halo.progenitor_next =construct_unique_id(tree_node->progenitor_next, i_tree);

  // Set group ids; needs to be modified for MPI
  tree_node->halo.group_halo_first=construct_unique_id(tree_node->group_halo_first,i_tree);
  tree_node->halo.group_halo_next =construct_unique_id(tree_node->group_halo_next, i_tree);
  
  // Walk the tree
  current=tree_node->progenitor_first;
  while(current!=NULL){
    assign_unique_ids_recursive(current,i_tree);
    current=current->progenitor_next;
  }
}

void compute_trees_vertical(char *filename_root_out,
                            char *filename_cat_root_in,
                            char *filename_snap_list_in,
                            int   n_files_groups,
                            int   n_files_subgroups,
                            int  *flag_clean){
  SID_fp      fp_in;
  SID_fp      fp_out;
  SID_fp      fp_out_MBP;
  char        filename_in[256];
  char        filename_out[256];
  char        filename_out_MBP[256];
  char        filename_output_file_root[256];
  char        filename_output_dir_horizontal[256];
  char        filename_output_dir_horizontal_trees[256];
  char        filename_output_dir_vertical[256];
  char        filename_output_dir_horizontal_groups_properties[256];
  char        filename_output_dir_horizontal_groups[256];
  char        filename_output_dir_horizontal_subgroups_properties[256];
  char        filename_output_dir_horizontal_subgroups[256];
  char        filename_output_vertical_root[256];
  char       *filename_output_dir_horizontal_properties;
  char       *line=NULL;
  size_t      line_length=0;
  int         i_bin;
  int         j_bin;
  int         n_snap;
  int         i_write;
  int         n_files;
  int         i_read;
  int         j_read;
  int         k_read;
  int         i_file;
  int         j_file;
  int         i_group;
  int         i_subgroup;
  int         j_subgroup;
  int         k_match;
  int         n_forests_local;
  int         i_rank;
  int         i_search;
  int         i_forest;
  int         i_tree;
  int         k_subgroup;
  int         group_descendant_id;
  int         subgroup_descendant_id;
  int         group_file_offset;
  int         group_file_index;
  int         forest_lo_rank;
  int         n_groups;
  int         n_subgroups;
  int         n_progenitors_max;
  int         n_trees_subgroup;
  int         n_trees_group;
  int         n_forests_subgroup;
  int         n_forests_group;
  int         n_tree_bins_groups;
  int         n_tree_bins_subgroups;
  int         n_tree_bins_groups_file;
  int         n_tree_bins_subgroups_file;
  int         n_tree_bins_groups_rank;
  int         n_tree_bins_subgroups_rank;
  int         n_subgroups_group;
  int         group_id;
  int         group_type;
  int         group_tree_id;
  int         subgroup_id;
  int         subgroup_type;
  int         subgroup_tree_id;
  int         min_sum;
  int         min_bin;
  int         n_write;
  int        *forest_lo_file=NULL;  
  int        *forest_hi_file=NULL;
  int        *tree_count_file=NULL;
  int         n_trees_file;
  int  flag_match_subgroups;
  int *i_forest_group;
  size_t *i_forest_group_index;
  int *i_forest_subgroup;
  size_t *i_forest_subgroup_index;
  int  i_forest_group_min;
  int  i_forest_subgroup_min;
  int *n_halos_tree_group;
  int *n_halos_tree_subgroup;
  int  n_halos_groups;
  int  n_halos_subgroups;
  int  n_halos_used;
  int  n_halos_target;
  int *n_halos_forest;
  int *forest_lo_group_file=NULL;  
  int *forest_hi_group_file=NULL;  
  int  forest_lo_group_local;  
  int  forest_hi_group_local;  
  int *forest_lo_group_rank=NULL;
  int *forest_hi_group_rank=NULL;
  int *tree_count_group_file=NULL; 
  int  tree_count_group_local;  
  int *tree_count_group_rank=NULL;
  int *forest_lo_subgroup_file=NULL;
  int *forest_hi_subgroup_file=NULL;
  int  forest_lo_subgroup_local;  
  int  forest_hi_subgroup_local;  
  int *forest_lo_subgroup_rank=NULL;
  int *forest_hi_subgroup_rank=NULL;
  int *tree_count_subgroup_file=NULL;
  int  tree_count_subgroup_local;
  int *tree_count_subgroup_rank=NULL;
  int  n_trees_subgroup_local;
  int  n_trees_group_local;
  int  subgroup_file_offset;
  int  subgroup_file_index;
  char group_text_prefix[4];
  int  n_conjoined;
  int  n_conjoined_total;
  halo_MBP_info     halo_MBP;
  tree_info       **trees;
  tree_node_info   *current=NULL;
  tree_node_info   *last   =NULL;
  tree_node_info   *next   =NULL;
  int               depth_first_index;
  int flag_write_init;
  int k_tree;
  int flag_init;
  int progenitor_score;
  int flag;
  int max_id=0;
  int n_halos_written;
  int halo_snap,descendant_snap;

  // Read the tree search/scan parameters
  int i_read_start;
  int i_read_stop;
  int i_read_step;
  int n_search;
  int flag_fix_bridges;
  int flag_compute_fragmented;
  int flag_compute_ghosts;
  read_tree_run_parameters(filename_root_out,
                           &i_read_start,
                           &i_read_stop,
                           &i_read_step,
                           &n_search,
                           &flag_fix_bridges,
                           &flag_compute_fragmented,
                           &flag_compute_ghosts);
  write_a_list(filename_snap_list_in,
               filename_root_out,
               i_read_start,
               i_read_stop,
               i_read_step);

  SID_log("Constructing vertical merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  // Assign trees to files
  SID_log("Assigning trees to files and ranks...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Compute a histogram of tree occupation by tree id
  for(i_read=i_read_stop,n_snap=0;i_read>=i_read_start;i_read-=i_read_step) n_snap++;

  // Determine the directory that the horizontal tree files are in
  strcpy(filename_output_file_root,filename_root_out);
  strip_path(filename_output_file_root);
  sprintf(filename_output_dir_horizontal,          "%s/horizontal",filename_root_out);
  sprintf(filename_output_dir_horizontal_groups,   "%s/groups",    filename_output_dir_horizontal);
  sprintf(filename_output_dir_horizontal_subgroups,"%s/subgroups", filename_output_dir_horizontal);
  sprintf(filename_output_dir_horizontal_trees,    "%s/trees",     filename_output_dir_horizontal);
  sprintf(filename_output_dir_vertical,            "%s/vertical",  filename_root_out);
  sprintf(filename_output_dir_horizontal_groups_properties,   "%s/properties",filename_output_dir_horizontal_groups);
  sprintf(filename_output_dir_horizontal_subgroups_properties,"%s/properties",filename_output_dir_horizontal_subgroups);
  mkdir(filename_output_dir_vertical,          02755);

  // Figure-out what the last-used snapshot is so we can read its header
  int i_read_last;
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step) i_read_last=i_read;

  // Determining iso-tree id mappings (onto tree structures, ranks and files).  Read header info
  //    from the i_read_start'th file since that is where the final tree tally gets written.
  int n_step_in;
  int n_search_in;
  int n_groups_max_in;
  int n_subgroups_max_in;
  sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read_last);
  SID_fopen(filename_in,"r",&fp_in);
  SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
  SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
  SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
  SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
  SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_in);
  SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_in);
  SID_fclose(&fp_in);
  n_progenitors_max=MAX(n_groups_max_in,n_subgroups_max_in);
  if(n_step_in!=i_read_step) SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);

  // Initialize arrays
  i_forest_group          =(int *)SID_malloc(sizeof(int)*n_trees_group);    // These are the linked tree ids that
  i_forest_subgroup       =(int *)SID_malloc(sizeof(int)*n_trees_subgroup); //   will be used for each iso-tree
  n_halos_tree_group      =(int *)SID_malloc(sizeof(int)*n_trees_group);
  n_halos_tree_subgroup   =(int *)SID_malloc(sizeof(int)*n_trees_subgroup);
  for(i_tree=0;i_tree<n_trees_group;i_tree++){
    i_forest_group[i_tree]    =i_tree;
    n_halos_tree_group[i_tree]=0;
  }
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
    i_forest_subgroup[i_tree]    =n_trees_subgroup+i_tree; // initialize to an impossibly high value
    n_halos_tree_subgroup[i_tree]=0;
  }

  // Determine the mapping between the horizontal tree IDs and the
  //   forest IDs that indicate how the isotrees are packaged into
  //   vertical trees.
  // Read-in all the tree halos and join trees into forests
  //    when substructure tree IDs start mixing between their initial values.
  //    Perform scan only over n_search to prevent over-linking.
  int n_trees_group_i;
  int n_trees_subgroup_i;
  int i_forest_group_i;
  int i_forest_subgroup_i;
  int n_halos_groups_unused=0;
  int n_halos_subgroups_unused=0;
  SID_log("Generating mapping of trees to forests...",SID_LOG_OPEN|SID_LOG_TIMER);
  for(i_read=i_read_stop,n_halos_groups=0,n_halos_subgroups=0;i_read>=i_read_start && j_read<n_search;i_read-=i_read_step){
    SID_log("Processing snapshot #%03d...",SID_LOG_OPEN,i_read);
    sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);
    SID_fopen(filename_in,"r",&fp_in);
    SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
    SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
    SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
    SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
    SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
    SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
    SID_fread_all(&n_trees_subgroup_i,sizeof(int),1,&fp_in);
    SID_fread_all(&n_trees_group_i,   sizeof(int),1,&fp_in);
    n_progenitors_max=MAX(n_groups_max_in,n_subgroups_max_in);
    if(n_step_in!=i_read_step) SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);
    for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
       SID_fread_all(&(group_id),           sizeof(int),1,&fp_in);
       SID_fread_all(&(group_type),         sizeof(int),1,&fp_in);
       SID_fread_all(&(group_descendant_id),sizeof(int),1,&fp_in);
       SID_fread_all(&(group_tree_id),      sizeof(int),1,&fp_in);
       SID_fread_all(&(group_file_offset),  sizeof(int),1,&fp_in);
       SID_fread_all(&(group_file_index),   sizeof(int),1,&fp_in);
       SID_fread_all(&(n_subgroups_group),  sizeof(int),1,&fp_in);
       if(group_id>=0 && group_tree_id>=0){
         i_forest_group_i=i_forest_group[group_tree_id];
         n_halos_tree_group[group_tree_id]++;
         n_halos_groups++;
       }
       else{
         i_forest_group_i=-1;
         n_halos_groups_unused++;
       }
       for(j_subgroup=0,i_forest_subgroup_min=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
         SID_fread_all(&(subgroup_id),           sizeof(int),1,&fp_in);
         SID_fread_all(&(subgroup_type),         sizeof(int),1,&fp_in);
         SID_fread_all(&(subgroup_descendant_id),sizeof(int),1,&fp_in);
         SID_fread_all(&(subgroup_tree_id),      sizeof(int),1,&fp_in);
         SID_fread_all(&(subgroup_file_offset),  sizeof(int),1,&fp_in);
         SID_fread_all(&(subgroup_file_index),   sizeof(int),1,&fp_in);
         if(subgroup_id>=0 && subgroup_tree_id>=0){
           i_forest_subgroup_i=i_forest_subgroup[subgroup_tree_id];
           n_halos_tree_subgroup[subgroup_tree_id]++;
           n_halos_subgroups++;
           // Initialize subgroup forest ID
           if(i_forest_subgroup_i>=n_trees_subgroup || i_forest_subgroup_i<0){
              i_forest_subgroup[subgroup_tree_id]=i_forest_group_i;
              i_forest_subgroup_i                =i_forest_subgroup[subgroup_tree_id];
           }
           // Join forests if this subgroup demands it.  Use the lower value of the two options.
           if(i_forest_subgroup_i<i_forest_group_i){
              i_forest_group[group_tree_id]=i_forest_subgroup_i;
              i_forest_group_i             =i_forest_group[group_tree_id];
           }
         }
         else
            n_halos_subgroups_unused++;
       }
    }
    SID_fclose(&fp_in);
    SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_log("Done.",SID_LOG_CLOSE);

  // Report some statistics
  SID_log("No. of    group isotrees = %d",SID_LOG_COMMENT,n_trees_group);
  SID_log("No. of subgroup isotrees = %d",SID_LOG_COMMENT,n_trees_subgroup);
  SID_log("No. of    groups         = %d",SID_LOG_COMMENT,n_halos_groups);
  SID_log("No. of subgroups         = %d",SID_LOG_COMMENT,n_halos_subgroups);
  SID_log("No. of    groups unused  = %d",SID_LOG_COMMENT,n_halos_groups_unused);
  SID_log("No. of subgroups unused  = %d",SID_LOG_COMMENT,n_halos_subgroups_unused);

  // Due to linking, there may be gaps in the forrest arrays at this point. Collapse them to fix this.

  // If an isotree has been joined to another to form a forest, one of them will
  //    have a reduced i_forest_group value (i_forest_subgroups point to this array at
  //    at this point, except for undefined subgroup trees which have 
  //    values > n_trees_group).  Follow the trail to the lowest value so that
  //    all isotrees meant for a single forrest have the same value.
  for(i_tree=0;i_tree<n_trees_group;i_tree++){
     int i_forest_old=i_forest_group[i_tree];
     if(i_forest_old>=0 && i_forest_old<n_trees_group)
        i_forest_group[i_tree]=i_forest_group[i_forest_old];
  }
  // At this point, the subgroups just point to group ids.  Give each
  //   subgroup to that group tree's forest.
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
     int i_forest_old=i_forest_subgroup[i_tree];
     if(i_forest_old>=0 && i_forest_old<n_trees_group)
        i_forest_subgroup[i_tree]=i_forest_group[i_forest_old];
  }
  
  //  ... eliminate subgroup gaps ...
  merge_sort(i_forest_subgroup,(size_t)n_trees_subgroup,&i_forest_subgroup_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  for(i_tree=0,k_tree=0;i_tree<n_trees_subgroup;){
    i_forest=i_forest_subgroup[i_forest_subgroup_index[i_tree]];
    while(i_forest_subgroup[i_forest_subgroup_index[i_tree]]==i_forest && i_tree<(n_trees_subgroup-1)) 
      i_forest_subgroup[i_forest_subgroup_index[i_tree++]]=k_tree;
    if(i_forest_subgroup[i_forest_subgroup_index[i_tree]]==i_forest) 
      i_forest_subgroup[i_forest_subgroup_index[i_tree++]]=k_tree;
    k_tree++;
  }
  n_forests_subgroup=k_tree;
  SID_free(SID_FARG i_forest_subgroup_index);

  //  ... then do groups ...
  int n_trees_forest_groups_max=0;
  merge_sort(i_forest_group,(size_t)n_trees_group,&i_forest_group_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  for(i_tree=0,k_tree=0;i_tree<n_trees_group;){
    int n_trees_forest_groups=0;
    i_forest=i_forest_group[i_forest_group_index[i_tree]];
    while(i_forest_group[i_forest_group_index[i_tree]]==i_forest && i_tree<(n_trees_group-1)){
      i_forest_group[i_forest_group_index[i_tree++]]=k_tree;
      n_trees_forest_groups++;
    }
    if(i_forest_group[i_forest_group_index[i_tree]]==i_forest){
      i_forest_group[i_forest_group_index[i_tree++]]=k_tree;
      n_trees_forest_groups++;
    }
    k_tree++;
    if(n_trees_forest_groups>n_trees_forest_groups_max)
       n_trees_forest_groups_max=n_trees_forest_groups;
  }
  n_forests_group=k_tree;
  SID_free(SID_FARG i_forest_group_index);

  // Turn the tree halo counts into forest halo counts
  int *n_halos_forest_group   =(int *)SID_calloc(sizeof(int)*n_forests_group);
  int *n_halos_forest_subgroup=(int *)SID_calloc(sizeof(int)*n_forests_subgroup);
  for(i_tree=0;i_tree<n_trees_group;i_tree++)
     n_halos_forest_group[i_forest_group[i_tree]]+=n_halos_tree_group[i_tree];
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++)
     n_halos_forest_subgroup[i_forest_group[i_tree]]+=n_halos_tree_subgroup[i_tree];
  SID_free(SID_FARG n_halos_tree_group);
  SID_free(SID_FARG n_halos_tree_subgroup);

  // Initialize some arrays
  forest_lo_group_rank    =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_hi_group_rank    =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  tree_count_group_rank   =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_lo_group_file    =(int *)SID_malloc(sizeof(int)*n_files_groups);
  forest_hi_group_file    =(int *)SID_malloc(sizeof(int)*n_files_groups);
  tree_count_group_file   =(int *)SID_malloc(sizeof(int)*n_files_groups);
  forest_lo_subgroup_rank =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_hi_subgroup_rank =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  tree_count_subgroup_rank=(int *)SID_malloc(sizeof(int)*SID.n_proc);  
  forest_lo_subgroup_file =(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
  forest_hi_subgroup_file =(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
  tree_count_subgroup_file=(int *)SID_malloc(sizeof(int)*n_files_subgroups);

  // Determine how subgroup trees will be distributed between ranks 
  for(i_rank=0,n_halos_used=0,i_tree=0;i_rank<SID.n_proc;i_rank++,i_tree++){
    forest_lo_subgroup_rank[i_rank]   =i_tree;
    forest_hi_subgroup_rank[i_rank]   =i_tree;
    tree_count_subgroup_rank[i_rank]=0;
    if(i_tree<n_forests_subgroup){
      n_halos_target                =(n_halos_subgroups-n_halos_used)/(SID.n_proc-i_rank);
      tree_count_subgroup_rank[i_rank]+=n_halos_forest_subgroup[i_tree];
      while(tree_count_subgroup_rank[i_rank]<n_halos_target && i_tree<(n_forests_subgroup-1)){
        i_tree++;
        forest_hi_subgroup_rank[i_rank]    =i_tree;
        tree_count_subgroup_rank[i_rank]+=n_halos_forest_subgroup[i_tree];
      }
      n_halos_used+=tree_count_subgroup_rank[i_rank];
    }
  }
  while(i_tree<(n_forests_subgroup-1)){
    i_tree++;
    forest_hi_subgroup_rank[SID.n_proc-1]    =i_tree;
    tree_count_subgroup_rank[SID.n_proc-1]+=n_halos_forest_subgroup[i_tree];
  }

  // Determine how group trees will be distributed between ranks 
  for(i_rank=0,n_halos_used=0,i_tree=0;i_rank<SID.n_proc;i_rank++,i_tree++){
    forest_lo_group_rank[i_rank]   =i_tree;
    forest_hi_group_rank[i_rank]   =i_tree;
    tree_count_group_rank[i_rank]=0;
    if(i_tree<n_forests_group){
      n_halos_target                =(n_halos_groups-n_halos_used)/(SID.n_proc-i_rank);
      tree_count_group_rank[i_rank]+=n_halos_forest_group[i_tree];
      while(tree_count_group_rank[i_rank]<n_halos_target && i_tree<(n_forests_group-1)){
        i_tree++;
        forest_hi_group_rank[i_rank]    =i_tree;
        tree_count_group_rank[i_rank]+=n_halos_forest_group[i_tree];
      }
      n_halos_used+=tree_count_group_rank[i_rank];
    }
  }
  while(i_tree<(n_forests_group-1)){
    i_tree++;
    forest_hi_group_rank[SID.n_proc-1]    =i_tree;
    tree_count_group_rank[SID.n_proc-1]+=n_halos_forest_group[i_tree];
  }

  // Determine how subgroup trees will be distributed between files
  for(i_file=0,n_halos_used=0,i_tree=0;i_file<n_files_subgroups;i_file++,i_tree++){
    forest_lo_subgroup_file[i_file]   =i_tree;
    forest_hi_subgroup_file[i_file]   =i_tree;
    tree_count_subgroup_file[i_file]=0;
    if(i_tree<n_forests_subgroup){
      n_halos_target                   =(n_halos_subgroups-n_halos_used)/(n_files_subgroups-i_file);
      tree_count_subgroup_file[i_file]+=n_halos_forest_subgroup[i_tree];
      while(tree_count_subgroup_file[i_file]<n_halos_target && i_tree<(n_forests_subgroup-1)){
        i_tree++;
        forest_hi_subgroup_file[i_file]    =i_tree;
        tree_count_subgroup_file[i_file]+=n_halos_forest_subgroup[i_tree];
      }
      n_halos_used+=tree_count_subgroup_file[i_file];
    }
  }
  while(i_tree<(n_forests_subgroup-1)){
    i_tree++;
    forest_hi_subgroup_file[n_files_subgroups-1]    =i_tree;
    tree_count_subgroup_file[n_files_subgroups-1]+=n_halos_forest_subgroup[i_tree];
  }

  // Determine how group trees will be distributed between files
  for(i_file=0,n_halos_used=0,i_tree=0;i_file<n_files_groups;i_file++,i_tree++){
    forest_lo_group_file[i_file]   =i_tree;
    forest_hi_group_file[i_file]   =i_tree;
    tree_count_group_file[i_file]=0;
    if(i_tree<n_forests_group){
      n_halos_target                =(n_halos_groups-n_halos_used)/(n_files_groups-i_file);
      tree_count_group_file[i_file]+=n_halos_forest_group[i_tree];
      while(tree_count_group_file[i_file]<n_halos_target && i_tree<(n_forests_group-1)){
        i_tree++;
        forest_hi_group_file[i_file]    =i_tree;
        tree_count_group_file[i_file]+=n_halos_forest_group[i_tree];
      }
      n_halos_used+=tree_count_group_file[i_file];
    }
  }
  while(i_tree<(n_forests_group-1)){
    i_tree++;
    forest_hi_group_file[n_files_groups-1]    =i_tree;
    tree_count_group_file[n_files_groups-1]+=n_halos_forest_group[i_tree];
  }

  // Set aliases for local values
  forest_lo_subgroup_local   =forest_lo_subgroup_rank[SID.My_rank];
  forest_hi_subgroup_local   =forest_hi_subgroup_rank[SID.My_rank];
  n_trees_subgroup_local   =forest_hi_subgroup_local-forest_lo_subgroup_local+1;
  tree_count_subgroup_local=tree_count_subgroup_rank[SID.My_rank];
  forest_lo_group_local      =forest_lo_group_rank[SID.My_rank];
  forest_hi_group_local      =forest_hi_group_rank[SID.My_rank];
  n_trees_group_local      =forest_hi_group_local-forest_lo_group_local+1;
  tree_count_group_local   =tree_count_group_rank[SID.My_rank];

  // Report decomposition results
  if(n_files_subgroups>1){
    for(i_bin=0;i_bin<n_files_subgroups;i_bin++)
      SID_log("Subgroup file #%4d will store %5d subgroup trees (%8d halos in total)",
               SID_LOG_COMMENT,i_bin+1,forest_hi_subgroup_file[i_bin]-forest_lo_subgroup_file[i_bin]+1,tree_count_subgroup_file[i_bin]);    
  }
  if(n_files_groups>1){
    for(i_bin=0;i_bin<n_files_groups;i_bin++)
      SID_log("Group    file #%4d will store %5d group    trees (%8d halos in total)",
              SID_LOG_COMMENT,i_bin+1,forest_hi_group_file[i_bin]-forest_lo_group_file[i_bin]+1,tree_count_group_file[i_bin]);    
  }
  if(SID.n_proc>1){
    for(i_bin=0;i_bin<SID.n_proc;i_bin++)
      SID_log("Rank     #%4d will process %5d subgroup trees (%8d halos in total)",
              SID_LOG_COMMENT,i_bin+1,forest_hi_subgroup_rank[i_bin]-forest_lo_subgroup_rank[i_bin]+1,tree_count_subgroup_rank[i_bin]);
    for(i_bin=0;i_bin<SID.n_proc;i_bin++)
      SID_log("Rank     #%4d will process %5d group    trees (%8d halos in total)",
              SID_LOG_COMMENT,i_bin+1,forest_hi_group_rank[i_bin]-forest_lo_group_rank[i_bin]+1,tree_count_group_rank[i_bin]);
  }
  SID_log("No. of    group forests  = %d",SID_LOG_COMMENT,n_forests_group);
  SID_log("No. of subgroup forests  = %d",SID_LOG_COMMENT,n_forests_subgroup);
  SID_log("Biggest forest size      = %d group isotrees",SID_LOG_COMMENT,n_trees_forest_groups_max);
  SID_log("Done.",SID_LOG_CLOSE);

  // VERTICAL TREE CONSTRUCTION STARTS HERE

  // Process subgroup trees (k_match==K_MATCH_SUBGROUPS) and then group trees (k_match==K_MATCH_GROUPS)
  //   If you want to change the order in which these are done, then reverse the preprocessor settings
  fp_catalog_info  fp_group_properties;
  fp_catalog_info  fp_subgroup_properties;
  halo_info       *group_properties;
  halo_info       *subgroup_properties;
  group_properties   =(halo_info *)SID_calloc(sizeof(halo_info));
  subgroup_properties=(halo_info *)SID_calloc(sizeof(halo_info));
  // Process subgroups and then groups
  for(k_match=0;k_match<2;k_match++){
    switch(k_match){
    case K_MATCH_GROUPS:
      filename_output_dir_horizontal_properties=filename_output_dir_horizontal_subgroups_properties;
      sprintf(filename_output_vertical_root,"%s",filename_output_dir_vertical);
      n_halos_forest   =n_halos_forest_subgroup;
      n_forests_local  =n_trees_subgroup_local;
      forest_lo_rank   =forest_lo_subgroup_rank[SID.My_rank];
      tree_count_file=tree_count_subgroup_file;
      forest_lo_file   =forest_lo_subgroup_file;
      forest_hi_file   =forest_hi_subgroup_file;
      n_write        =n_files_subgroups;
      sprintf(group_text_prefix,"sub");
      break;
    case K_MATCH_SUBGROUPS:
      filename_output_dir_horizontal_properties=filename_output_dir_horizontal_groups_properties;
      sprintf(filename_output_vertical_root,"%s",filename_output_dir_vertical);
      n_halos_forest   =n_halos_forest_group;
      n_forests_local  =n_trees_group_local;
      forest_lo_rank   =forest_lo_group_rank[SID.My_rank];
      tree_count_file=tree_count_group_file;
      forest_lo_file   =forest_lo_group_file;
      forest_hi_file   =forest_hi_group_file;
      n_write        =n_files_groups;
      sprintf(group_text_prefix,"");
      break;      
    }
    
    // Skip this iteration if the number of
    //   output files is set to <=0
    if(n_write<=0)
      continue;

    // Initialize tree arrays
    trees=(tree_info **)SID_malloc(sizeof(tree_info *)*n_forests_local);
    for(i_tree=0;i_tree<n_forests_local;i_tree++)
      init_tree(n_snap,&trees[i_tree]);
    
    // Print a status message
    SID_log("Rendering %sgroup trees vertical...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
    if((*flag_clean))
      SID_log("(horizontal tree files will be removed)...",SID_LOG_CONTINUE);

    // Loop over all the horizontal tree files in order of decreasing snapshot number, hanging halos on the trees as we go
    int n_subgroups_orphaned=0;
    for(i_read=i_read_stop,i_file=n_snap-1,flag_init=TRUE;i_read>=i_read_start;i_read-=i_read_step,i_file--,flag_init=FALSE){
      SID_log("Processing snapshot %03d (%03d of %03d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,i_file+1,n_snap);

      // This counter makes sure that the halo snap index in the trees
      //   is continuous, even if we are skipping snapshots
      halo_snap=i_file;

      // Initialize this snapshot's halo list here
      for(i_tree=0;i_tree<n_forests_local;i_tree++){
        trees[i_tree]->neighbour_halos[i_file]=NULL;
        trees[i_tree]->n_neighbours[i_file]   =0;
      }
      sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);

      // Open properties catalog
      fopen_catalog(filename_cat_root_in,
                    i_read,
                    READ_CATALOG_GROUPS|READ_CATALOG_PROPERTIES,
                    &fp_group_properties);
      if(k_match==K_MATCH_SUBGROUPS)
         fopen_catalog(filename_cat_root_in,
                       i_read,
                       READ_CATALOG_SUBGROUPS|READ_CATALOG_PROPERTIES,
                       &fp_subgroup_properties);

      // Open horizontal tree file
      SID_fopen(filename_in,"r",&fp_in);

      // Read header
      SID_fread_all(&n_step_in,         sizeof(int),1,&fp_in);
      SID_fread_all(&n_search_in,       sizeof(int),1,&fp_in);
      SID_fread_all(&n_groups,          sizeof(int),1,&fp_in);
      SID_fread_all(&n_subgroups,       sizeof(int),1,&fp_in);
      SID_fread_all(&n_groups_max_in,   sizeof(int),1,&fp_in);
      SID_fread_all(&n_subgroups_max_in,sizeof(int),1,&fp_in);
      SID_fread_all(&n_trees_subgroup,  sizeof(int),1,&fp_in);
      SID_fread_all(&n_trees_group,     sizeof(int),1,&fp_in);
      n_progenitors_max=MAX(n_groups_max_in,n_subgroups_max_in);
      if(n_step_in!=i_read_step) 
         SID_trap_error("Snapshot step sizes don't match (ie. %d!=%d)",ERROR_LOGIC,n_step_in,i_read_step);
      
      // Read each group in turn
      for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){

        // Read horizontal trees for groups
        SID_fread_all(&(group_id),           sizeof(int),1,&fp_in);
        SID_fread_all(&(group_type),         sizeof(int),1,&fp_in);
        SID_fread_all(&(group_descendant_id),sizeof(int),1,&fp_in);
        SID_fread_all(&(group_tree_id),      sizeof(int),1,&fp_in);
        SID_fread_all(&(group_file_offset),  sizeof(int),1,&fp_in);
        SID_fread_all(&(group_file_index),   sizeof(int),1,&fp_in);
        SID_fread_all(&(n_subgroups_group),  sizeof(int),1,&fp_in);
        if(group_tree_id>=0)
          group_tree_id=i_forest_group[group_tree_id];
        else
          group_tree_id=-1;

        // Read halo information from catalog files (needed even for subgroup trees; we need FoF masses for the most massive progenitors)
        fread_catalog_file(&fp_group_properties,group_properties,NULL,NULL,i_group);

        // If we are processing subgroup trees ...
        if(k_match==K_MATCH_SUBGROUPS){
          // Read each subgroup in turn
          for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
            SID_fread_all(&(subgroup_id),           sizeof(int),1,&fp_in);
            SID_fread_all(&(subgroup_type),         sizeof(int),1,&fp_in);
            SID_fread_all(&(subgroup_descendant_id),sizeof(int),1,&fp_in);
            SID_fread_all(&(subgroup_tree_id),      sizeof(int),1,&fp_in);
            SID_fread_all(&(subgroup_file_offset),  sizeof(int),1,&fp_in);
            SID_fread_all(&(subgroup_file_index),   sizeof(int),1,&fp_in);
            // Ignore negative ids
            if(subgroup_tree_id>=0)
              subgroup_tree_id=i_forest_subgroup[subgroup_tree_id];
            else
              subgroup_tree_id=-1;
            if(subgroup_id>=0 && subgroup_tree_id>=0){
              // If this subgroup belongs to a local tree ...
              i_tree=subgroup_tree_id-forest_lo_rank;
              if(i_tree>=0 && i_tree<n_forests_local){ 
                // ... create a new branch and add it to its tree ...
                if(subgroup_file_offset<=0)
                  descendant_snap=-1; // Needed for halos in the root snapshot
                else
                  descendant_snap=(i_file+subgroup_file_offset);
                // ... read halo information from catalog files ...
                fread_catalog_file(&fp_subgroup_properties,subgroup_properties,NULL,NULL,i_subgroup);
                // ... set the most massive progenitor's mass to the FoF mass ...
                if(j_subgroup==0)
                   subgroup_properties->M_vir=group_properties->M_vir;
                // ... adjust snap counter to make things work with skipped snaps ...
                subgroup_properties->snap_num=halo_snap;
                // ... add this halo to the trees ...
                int flag_found_group=
                add_node_to_tree(trees[i_tree],
                                 subgroup_type,
                                 subgroup_id,
                                 group_id,
                                 subgroup_descendant_id,
                                 halo_snap,
                                 descendant_snap,
                                 subgroup_properties);
                if(j_subgroup>0 && !flag_found_group)
                   n_subgroups_orphaned++;
              }
            }
          }
        }
        // ... else, process group trees 
        else{
          SID_fseek(&fp_in,sizeof(int),n_subgroups_group*6,SID_SEEK_CUR);
          // Ignore negative IDs
          if(group_id>=0 && group_tree_id>=0){
            i_tree=group_tree_id-forest_lo_rank;
            // If this group belongs to a local tree ...
            if(i_tree>=0 && i_tree<n_forests_local){ 
              // ... if so, create a new branch and add it to the tree ...
              if(group_file_offset<=0)
                descendant_snap=-1; // Needed for halos in the root snapshot
              else
                descendant_snap=(i_file+group_file_offset);
              // ... adjust snap counter to make things work with skipped snaps ...
              group_properties->snap_num=halo_snap;
              // ... add this halo to the trees ...
              add_node_to_tree(trees[i_tree],
                               group_type,
                               group_id,
                               group_id,
                               group_descendant_id,
                               halo_snap,
                               descendant_snap,                      
                               group_properties);
            }
          }
        } // if k_match
      } // i_group
      SID_fclose(&fp_in);
      fclose_catalog(&fp_group_properties);
      if(k_match==K_MATCH_SUBGROUPS)
         fclose_catalog(&fp_subgroup_properties);

      // If flag_clean=TRUE, then delete the input files used here.
      //if((*flag_clean)==TRUE){
      //  if(k_match==1){
      //    sprintf(filename_in,"%s/horizontal_trees_%03d.dat",filename_output_dir_horizontal_trees,i_read);
      //    remove(filename_in);
      //  }
      //}      
      SID_log("Done.",SID_LOG_CLOSE);
    } // i_read
    SID_log("Done.",SID_LOG_CLOSE);

    // Report any subgroups that become orphaned from their parent halo
    if(k_match==K_MATCH_SUBGROUPS){
       SID_Allreduce(SID_IN_PLACE,&n_subgroups_orphaned,1,SID_INT,SID_SUM,SID.COMM_WORLD);
       SID_log("n_subgroups_orphaned= %d",SID_LOG_COMMENT,n_subgroups_orphaned);
    }

    // Finalize trees
    finalize_trees_vertical(trees,&(n_halos_forest[forest_lo_rank]),n_forests_local,n_snap,TREE_PROGENITOR_ORDER_DELUCIA);

    // Write trees
    write_trees_vertical(trees,&(n_halos_forest[forest_lo_rank]),n_forests_local,forest_lo_file,forest_hi_file,tree_count_file,n_write,filename_output_vertical_root,group_text_prefix);

    // Free trees
    for(i_tree=0;i_tree<n_forests_local;i_tree++)
      free_tree(&trees[i_tree]);
    SID_free((void **)&trees);

  } // k_match
  SID_free(SID_FARG group_properties);
  SID_free(SID_FARG subgroup_properties);
  
  // Clean-up
  SID_free(SID_FARG line);
  SID_free(SID_FARG i_forest_group);
  SID_free(SID_FARG i_forest_subgroup);
  SID_free(SID_FARG n_halos_forest_group);
  SID_free(SID_FARG n_halos_forest_subgroup);
  SID_free(SID_FARG forest_lo_group_rank);
  SID_free(SID_FARG forest_hi_group_rank);
  SID_free(SID_FARG tree_count_group_rank);
  SID_free(SID_FARG forest_lo_subgroup_rank);
  SID_free(SID_FARG forest_hi_subgroup_rank);
  SID_free(SID_FARG tree_count_subgroup_rank);
  SID_free(SID_FARG forest_lo_group_file);
  SID_free(SID_FARG forest_hi_group_file);
  SID_free(SID_FARG tree_count_group_file);
  SID_free(SID_FARG forest_lo_subgroup_file);
  SID_free(SID_FARG forest_hi_subgroup_file);
  SID_free(SID_FARG tree_count_subgroup_file);

  SID_log("Done.",SID_LOG_CLOSE);
}

