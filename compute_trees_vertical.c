#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void init_tree(int n_search,tree_info **tree);
void init_tree(int n_search,tree_info **tree){
  int i_search;
  (*tree)                     =(tree_info       *)SID_malloc(sizeof(tree_info));
  (*tree)->n_neighbours       =(int             *)SID_malloc(sizeof(int)*n_search);
  (*tree)->neighbour_halos    =(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_search);
  (*tree)->neighbour_halo_last=(tree_node_info **)SID_malloc(sizeof(tree_node_info *)*n_search);
  for(i_search=0;i_search<n_search;i_search++){
    (*tree)->n_neighbours[i_search]       =0;
    (*tree)->neighbour_halos[i_search]    =NULL;
    (*tree)->neighbour_halo_last[i_search]=NULL;
  }
  (*tree)->root     =NULL;
  (*tree)->last_leaf=NULL;
}

void free_tree(tree_info **tree);
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

void add_node_to_tree(tree_info  *tree,
                      int         halo_id,
                      int         group_id,
                      int         descendant_id,
                      int         halo_snap,
                      int         descendant_snap,                      
                      halo_info  *properties);
void add_node_to_tree(tree_info  *tree,
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
  new_node->depth_first_index    =       -1; // Default; back-filled later
  new_node->n_progenitors        =        0;
  new_node->group_id             = group_id;
  new_node->halo_id              =  halo_id;
  new_node->descendant           =     NULL; // Set below
  new_node->progenitor_first     =     NULL; // Set below
  new_node->progenitor_next      =     NULL; // Set below
  new_node->progenitor_last      =     NULL; // Set below
  new_node->group_halo_first     = new_node; // Set below if group_id>=0
  new_node->group_halo_next      =     NULL; // Set below
  new_node->neighbour_halo_next  =     NULL; // Set below
  new_node->next                 =     NULL;

  // Find the halo's descendant and set various pointers
  //   All used halos must be added to the neighbour list (done below) for this to work
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
        while(next_node!=NULL){
          last_node=next_node;
          next_node=last_node->group_halo_next;
        }
        last_node->group_halo_next=new_node;        
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
}

int construct_unique_id(tree_node_info *tree_node,int tree_number);
int construct_unique_id(tree_node_info *tree_node,int tree_number){
  if(tree_node!=NULL){
    //return(tree_number*1000000+tree_node->depth_first_index);
    return(tree_node->depth_first_index);
  }
  else
    return(-1);    
}

void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode);
void assign_progenitor_order_recursive(tree_node_info *tree,int *M_i,int mode){
  tree_node_info  *first_new;
  tree_node_info  *last_new;
  tree_node_info  *current;
  tree_node_info **progenitors;
  int              i_progenitor;
  size_t          *M_iN_index;
  size_t          *M_iN_index_index; 
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
    merge_sort(M_iN,      (size_t)tree->n_progenitors,&M_iN_index,      SID_INT,   SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    merge_sort(M_iN_index,(size_t)tree->n_progenitors,&M_iN_index_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    first_new=progenitors[M_iN_index[tree->n_progenitors-1]]; 
    last_new =progenitors[M_iN_index[0]]; 
    tree->progenitor_first=first_new;
    tree->progenitor_last =last_new;
    for(i_progenitor=0;i_progenitor<tree->n_progenitors;i_progenitor++){
      if(progenitors[i_progenitor]!=last_new)
        progenitors[i_progenitor]->progenitor_next=progenitors[M_iN_index[M_iN_index_index[i_progenitor]-1]]; 
      else
        progenitors[i_progenitor]->progenitor_next=NULL;
    }

    /*
    if(tree->n_progenitors>1){
      i_progenitor=0;
      current     =tree->progenitor_first;
      while(current!=NULL){
        fprintf(stderr,"i=%4d pf=%p pl=%p f=%p c=%p n=%p n_p=%d M_vir=%le\n",i_progenitor++,tree->progenitor_first,tree->progenitor_last,current->progenitor_last,current,current->progenitor_next,current->halo.n_particles,current->halo.M_vir);
        current=current->progenitor_next;
      }
    }
    */

    SID_free(SID_FARG progenitors);
    SID_free(SID_FARG M_iN);
    SID_free(SID_FARG M_iN_index);
    SID_free(SID_FARG M_iN_index_index);
  }

  // Add this progenitor's score to the descendant's sum (see De Lucia and Blaizot (2006))
  (*M_i)=N_i+max_M_iN;
}

void assign_group_halo_order(tree_info *tree,int i_snap);
void assign_group_halo_order(tree_info *tree,int i_snap){
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
  int             *halo_size;
  size_t          *group_ids_index;
  size_t          *halo_size_index;
  size_t          *halo_size_index_index;

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
        SID_trap_error("There's a problem with the number of neighbours in assign_group_halo_order (n_neighbours=%d)!",ERROR_LOGIC,n_neighbours);
      neighbours[i_halo]=current;
      group_ids[i_halo] =current->group_id;
      current           =current->neighbour_halo_next;
      i_halo++;
    }
    merge_sort(group_ids,(size_t)n_neighbours,&group_ids_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

    // Find size of largest group and initialize a temporary array
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

    // Scan neighbour list (now ordered by ID) and order by size within each group
    if(largest_group>1){
      halo_size=(int *)SID_malloc(sizeof(int)*largest_group);
      i_halo   =0;
      while(i_halo<n_neighbours){
        // Create a list of halo sizes for each group
        n_in_group=0;
        while(group_ids[group_ids_index[i_halo+n_in_group]]==group_ids[group_ids_index[i_halo]]){
          halo_size[n_in_group]=neighbours[group_ids_index[i_halo+n_in_group]]->halo.n_particles;
          n_in_group++;
          if((i_halo+n_in_group)>=n_neighbours)
            break;
        }
        // Correct the ordering within each group ...
        if(n_in_group>1){
          // ... sort halo sizes (ascending) ...
          merge_sort(halo_size,      (size_t)n_in_group,&halo_size_index,      SID_INT,   SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          merge_sort(halo_size_index,(size_t)n_in_group,&halo_size_index_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          // ... set the new pointers ...
          new_first=neighbours[group_ids_index[i_halo+halo_size_index[n_in_group-1]]];
          new_last =neighbours[group_ids_index[i_halo+halo_size_index[0]]];
          for(j_halo=i_halo,k_halo=0;k_halo<n_in_group;j_halo++,k_halo++){
            if(neighbours[group_ids_index[j_halo]]!=new_last)
              neighbours[group_ids_index[j_halo]]->group_halo_next=neighbours[group_ids_index[i_halo+halo_size_index[halo_size_index_index[k_halo]-1]]];
            else
              neighbours[group_ids_index[j_halo]]->group_halo_next=NULL;
            neighbours[group_ids_index[j_halo]]->group_halo_first=new_first;
          }
          SID_free(SID_FARG halo_size_index);
          SID_free(SID_FARG halo_size_index_index);
        }
        i_halo+=n_in_group;
      }
      SID_free(SID_FARG halo_size);
    }

    SID_free(SID_FARG group_ids);
    SID_free(SID_FARG group_ids_index);
    SID_free(SID_FARG neighbours);
  }
}

void assign_depth_first_index_recursive(tree_node_info *tree,int *depth_first_index);
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

void assign_unique_ids_recursive(tree_node_info *tree_node,int i_tree);
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

int write_tree_vertical_halos_recursive(tree_node_info *tree_node,SID_fp *fp_out,SID_fp *fp_out_MBP);
int write_tree_vertical_halos_recursive(tree_node_info *tree_node,SID_fp *fp_out,SID_fp *fp_out_MBP){
  tree_node_info *current;
  halo_info      *halo;
  halo_MBP_info   halo_MBP;
  int             n_halos_written=0;

  // Write tree halos
  halo=&(tree_node->halo);
  SID_fwrite(halo,sizeof(halo_info),1,fp_out);

  // Write MBPs
  if(fp_out_MBP!=NULL){
    halo_MBP.most_bound_id   =halo->most_bound_id;
    halo_MBP.snap_num        =halo->snap_num;
    halo_MBP.halo_index      =halo->halo_index;
    halo_MBP.group_halo_first=halo->group_halo_first;
    halo_MBP.pos[0]          =halo->pos[0];
    halo_MBP.pos[1]          =halo->pos[1];
    halo_MBP.pos[2]          =halo->pos[2];
    halo_MBP.vel[0]          =halo->vel[0];
    halo_MBP.vel[1]          =halo->vel[1];
    halo_MBP.vel[2]          =halo->vel[2];
    SID_fwrite(&halo_MBP,sizeof(halo_MBP_info),1,fp_out_MBP);
  }
  
  n_halos_written++;
  current=tree_node->progenitor_first;
  while(current!=NULL){
    n_halos_written+=write_tree_vertical_halos_recursive(current,fp_out,fp_out_MBP);
    current=current->progenitor_next;
  }
  return(n_halos_written);
}

void compute_trees_vertical(char *filename_root_out,
                            int   i_read_start,
                            int   i_read_stop,
                            int   i_read_step,
                            int   n_search,
                            int   n_files_groups,
                            int   n_files_subgroups,
                            int  *flag_clean){
  FILE       *fp_in;
  SID_fp      fp_out;
  SID_fp      fp_out_MBP;
  SID_fp      fp_properties;
  char        filename_in[256];
  char        filename_out[256];
  char        filename_out_MBP[256];
  char        filename_properties[256];
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
  int         n_trees_local;
  int         i_rank;
  int         i_search;
  int         i_tree;
  int         j_tree;
  int         k_subgroup;
  int         group_descendant_id;
  int         subgroup_descendant_id;
  int         group_file_offset;
  int         tree_lo_rank;
  int         n_groups;
  int         n_subgroups;
  int         n_progenitors_max;
  int         n_trees_subgroup;
  int         n_trees_group;
  int         n_tree_bins_groups;
  int         n_tree_bins_subgroups;
  int         n_tree_bins_groups_file;
  int         n_tree_bins_subgroups_file;
  int         n_tree_bins_groups_rank;
  int         n_tree_bins_subgroups_rank;
  int         n_subgroups_group;
  int         group_id;
  int         group_tree_id;
  int         subgroup_id;
  int         subgroup_tree_id;
  int         min_sum;
  int         test_sum;
  int         min_bin;
  int         n_write;
  int        *tree_lo_file=NULL;  
  int        *tree_hi_file=NULL;
  int        *tree_count_file=NULL;
  int         n_trees_file;
  int  flag_match_subgroups;
  int *i_tree_group;
  size_t *i_tree_group_index;
  int *i_tree_subgroup;
  size_t *i_tree_subgroup_index;
  int  i_tree_group_min;
  int  i_tree_subgroup_min;
  int *n_halos_tree_group;
  int *n_halos_tree_subgroup;
  int  n_halos_groups;
  int  n_halos_subgroups;
  int  n_halos_used;
  int  n_halos_target;
  int *n_halos_tree;
  int *tree_lo_group_file=NULL;  
  int *tree_hi_group_file=NULL;  
  int  tree_lo_group_local;  
  int  tree_hi_group_local;  
  int *tree_lo_group_rank=NULL;
  int *tree_hi_group_rank=NULL;
  int *tree_count_group_file=NULL; 
  int  tree_count_group_local;  
  int *tree_count_group_rank=NULL;
  int *tree_lo_subgroup_file=NULL;
  int *tree_hi_subgroup_file=NULL;
  int  tree_lo_subgroup_local;  
  int  tree_hi_subgroup_local;  
  int *tree_lo_subgroup_rank=NULL;
  int *tree_hi_subgroup_rank=NULL;
  int *tree_count_subgroup_file=NULL;
  int  tree_count_subgroup_local;
  int *tree_count_subgroup_rank=NULL;
  int  n_trees_subgroup_local;
  int  n_trees_group_local;
  int  subgroup_file_offset;
  char group_text_prefix[4];
  int  n_conjoined;
  int  n_conjoined_total;
  halo_info         properties;
  halo_info        *halo;
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

  SID_log("Constructing vertical merger trees for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  // Assign trees to files
  SID_log("Assigning trees to files and ranks...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Compute a histogram of tree occupation by tree id
  for(i_read=i_read_stop,n_snap=0;i_read>=i_read_start;i_read-=i_read_step) n_snap++;

  n_search=n_search+2;
  
  // Determining iso-tree id mappings (onto tree structures, ranks and files)
  //   Loop over each snapshot ...
  for(i_read=i_read_stop,flag_init=TRUE,n_halos_groups=0,n_halos_subgroups=0;i_read>=i_read_start;i_read-=i_read_step){
    sprintf(filename_in,"%s.trees_horizontal_%d",filename_root_out,i_read);
    fp_in=fopen(filename_in,"r");
    grab_next_line(fp_in,&line,&line_length);
    grab_int(line,1,&n_groups);
    grab_int(line,2,&n_subgroups);
    grab_int(line,4,&n_trees_subgroup);
    grab_int(line,5,&n_trees_group);
    // Initialize the mapping of group and subgroup halo tree identities onto the final tree structure identities
    if(flag_init){
      // Allocate arrays
      i_tree_group            =(int *)SID_malloc(sizeof(int)*n_trees_group);    // These are the linked tree ids that
      i_tree_subgroup         =(int *)SID_malloc(sizeof(int)*n_trees_subgroup); //   will be used for each iso-tree
      n_halos_tree_group      =(int *)SID_malloc(sizeof(int)*n_trees_group);
      n_halos_tree_subgroup   =(int *)SID_malloc(sizeof(int)*n_trees_subgroup);
      tree_lo_group_rank      =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
      tree_hi_group_rank      =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
      tree_count_group_rank   =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
      tree_lo_group_file      =(int *)SID_malloc(sizeof(int)*n_files_groups);
      tree_hi_group_file      =(int *)SID_malloc(sizeof(int)*n_files_groups);
      tree_count_group_file   =(int *)SID_malloc(sizeof(int)*n_files_groups);
      tree_lo_subgroup_rank   =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
      tree_hi_subgroup_rank   =(int *)SID_malloc(sizeof(int)*SID.n_proc);  
      tree_count_subgroup_rank=(int *)SID_malloc(sizeof(int)*SID.n_proc);  
      tree_lo_subgroup_file   =(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
      tree_hi_subgroup_file   =(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
      tree_count_subgroup_file=(int *)SID_malloc(sizeof(int)*n_files_subgroups);  
      // Initialize arrays
      for(i_tree=0;i_tree<n_trees_group;i_tree++){
        i_tree_group[i_tree]      =n_trees_group+1;
        n_halos_tree_group[i_tree]=0;
      }
      for(i_tree=0;i_tree<n_trees_subgroup;i_tree++){
        i_tree_subgroup[i_tree]      =n_trees_subgroup+1;
        n_halos_tree_subgroup[i_tree]=0;
      }
      // Process the first file; initialize tree ids to final group indices
      for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
        grab_next_line(fp_in,&line,&line_length);
        grab_int(line,1,&group_id);
        grab_int(line,3,&group_tree_id);
        grab_int(line,5,&n_subgroups_group);
        if(group_id>=0 && group_tree_id>=0){
          i_tree_group[group_tree_id]=group_tree_id;
          n_halos_groups++;
        }
        for(i_subgroup=0;i_subgroup<n_subgroups_group;i_subgroup++){
          grab_int(line,6+i_subgroup*4,&subgroup_id);
          grab_int(line,8+i_subgroup*4,&subgroup_tree_id);
          if(subgroup_id>=0 && subgroup_tree_id>=0){
            if(group_tree_id>=0)
              i_tree_subgroup[subgroup_tree_id]=group_tree_id;
            n_halos_subgroups++;
          }
        }
      }
      flag_init=FALSE;
    }
    // Loop over all the other snapshots, linking interwoven trees as we go
    else{
      for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
        grab_next_line(fp_in,&line,&line_length);
        grab_int(line,1,&group_id);
        grab_int(line,3,&group_tree_id);
        grab_int(line,5,&n_subgroups_group);
        if(group_id>=0 && group_tree_id>=0)
          n_halos_groups++;
        for(i_subgroup=0;i_subgroup<n_subgroups_group;i_subgroup++){
          grab_int(line,6+i_subgroup*4,&subgroup_id);
          grab_int(line,8+i_subgroup*4,&subgroup_tree_id);
          if(subgroup_id>=0 && subgroup_tree_id>=0){
            if(group_tree_id>=0)
              i_tree_subgroup[subgroup_tree_id]=group_tree_id;
            n_halos_subgroups++;
          }
        }
      }
    }
    fclose(fp_in);
  }
  SID_log("(n_groups=%d n_subgroups=%d)...",SID_LOG_CONTINUE,n_halos_groups,n_halos_subgroups);
  
  // Due to linking, there may be gaps in the i_tree arrays at this point. Compress them to fix this.
  merge_sort(i_tree_subgroup,(size_t)n_trees_subgroup,&i_tree_subgroup_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  for(i_tree=0,k_tree=0;i_tree<n_trees_subgroup;){
    j_tree=i_tree_subgroup[i_tree_subgroup_index[i_tree]];
    while(i_tree_subgroup[i_tree_subgroup_index[i_tree]]==j_tree && i_tree<n_trees_subgroup-1) 
      i_tree_subgroup[i_tree_subgroup_index[i_tree++]]=k_tree;
    if(i_tree_subgroup[i_tree_subgroup_index[i_tree]]==j_tree) 
      i_tree_subgroup[i_tree_subgroup_index[i_tree++]]=k_tree;
    k_tree++;
  }
  n_trees_subgroup=k_tree;
  SID_free((void **)&i_tree_subgroup_index);
  merge_sort(i_tree_group,(size_t)n_trees_group,&i_tree_group_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
  for(i_tree=0,k_tree=0;i_tree<n_trees_group;){
    j_tree=i_tree_group[i_tree_group_index[i_tree]];
    while(i_tree_group[i_tree_group_index[i_tree]]==j_tree && i_tree<n_trees_group-1) 
      i_tree_group[i_tree_group_index[i_tree++]]=k_tree;
    if(i_tree_group[i_tree_group_index[i_tree]]==j_tree) 
      i_tree_group[i_tree_group_index[i_tree++]]=k_tree;
    k_tree++;
  }
  n_trees_group=k_tree;
  SID_free((void **)&i_tree_group_index);

  // Determine tree halo counts
  for(i_read=i_read_stop;i_read>=i_read_start;i_read-=i_read_step){
    sprintf(filename_in,"%s.trees_horizontal_%d",filename_root_out,i_read);
    fp_in=fopen(filename_in,"r");
    grab_next_line(fp_in,&line,&line_length);
    grab_int(line,1,&n_groups);
    for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
      grab_next_line(fp_in,&line,&line_length);
      grab_int(line,1,&group_id);
      grab_int(line,3,&group_tree_id);
      grab_int(line,5,&n_subgroups_group);
      if(group_id>=0 && group_tree_id>=0)
        n_halos_tree_group[i_tree_group[group_tree_id]]++;
      for(i_subgroup=0,i_tree_subgroup_min=0;i_subgroup<n_subgroups_group;i_subgroup++){
        grab_int(line,6+i_subgroup*4,&subgroup_id);
        grab_int(line,8+i_subgroup*4,&subgroup_tree_id);
        if(subgroup_id>=0 && subgroup_tree_id>=0)
          n_halos_tree_subgroup[i_tree_subgroup[subgroup_tree_id]]++;
      }
    }
    fclose(fp_in);
  }

  // Determine the number of trees we are left with now
  for(i_tree=0;i_tree<n_trees_subgroup;i_tree++)
    if(n_halos_tree_subgroup[i_tree]>0) j_tree=i_tree;
  n_trees_subgroup=j_tree+1;
  for(i_tree=0;i_tree<n_trees_group;i_tree++)
    if(n_halos_tree_group[i_tree]>0) j_tree=i_tree;
  n_trees_group=j_tree+1;

  // Determine how subgroup trees will be distributed between ranks 
  for(i_rank=0,n_halos_used=0,i_tree=0;i_rank<SID.n_proc;i_rank++,i_tree++){
    tree_lo_subgroup_rank[i_rank]   =i_tree;
    tree_hi_subgroup_rank[i_rank]   =i_tree;
    tree_count_subgroup_rank[i_rank]=0;
    if(i_tree<n_trees_subgroup){
      n_halos_target                =(n_halos_subgroups-n_halos_used)/(SID.n_proc-i_rank);
      tree_count_subgroup_rank[i_rank]+=n_halos_tree_subgroup[i_tree];
      while(tree_count_subgroup_rank[i_rank]<n_halos_target && i_tree<n_trees_subgroup-1){
        i_tree++;
        tree_hi_subgroup_rank[i_rank]    =i_tree;
        tree_count_subgroup_rank[i_rank]+=n_halos_tree_subgroup[i_tree];
      }
      n_halos_used+=tree_count_subgroup_rank[i_rank];
    }
  }
  while(i_tree<n_trees_subgroup-1){
    i_tree++;
    tree_hi_subgroup_rank[SID.n_proc-1]    =i_tree;
    tree_count_subgroup_rank[SID.n_proc-1]+=n_halos_tree_subgroup[i_tree];
  }

  // Determine how group trees will be distributed between ranks 
  for(i_rank=0,n_halos_used=0,i_tree=0;i_rank<SID.n_proc;i_rank++,i_tree++){
    tree_lo_group_rank[i_rank]   =i_tree;
    tree_hi_group_rank[i_rank]   =i_tree;
    tree_count_group_rank[i_rank]=0;
    if(i_tree<n_trees_group){
      n_halos_target                =(n_halos_groups-n_halos_used)/(SID.n_proc-i_rank);
      tree_count_group_rank[i_rank]+=n_halos_tree_group[i_tree];
      while(tree_count_group_rank[i_rank]<n_halos_target && i_tree<n_trees_group-1){
        i_tree++;
        tree_hi_group_rank[i_rank]    =i_tree;
        tree_count_group_rank[i_rank]+=n_halos_tree_group[i_tree];
      }
      n_halos_used+=tree_count_group_rank[i_rank];
    }
  }
  while(i_tree<n_trees_group-1){
    i_tree++;
    tree_hi_group_rank[SID.n_proc-1]    =i_tree;
    tree_count_group_rank[SID.n_proc-1]+=n_halos_tree_group[i_tree];
  }

  // Determine how subgroup trees will be distributed between files
  for(i_file=0,n_halos_used=0,i_tree=0;i_file<n_files_subgroups;i_file++,i_tree++){
    tree_lo_subgroup_file[i_file]   =i_tree;
    tree_hi_subgroup_file[i_file]   =i_tree;
    tree_count_subgroup_file[i_file]=0;
    if(i_tree<n_trees_subgroup){
      n_halos_target                   =(n_halos_subgroups-n_halos_used)/(n_files_subgroups-i_file);
      tree_count_subgroup_file[i_file]+=n_halos_tree_subgroup[i_tree];
      while(tree_count_subgroup_file[i_file]<n_halos_target && i_tree<n_trees_subgroup-1){
        i_tree++;
        tree_hi_subgroup_file[i_file]    =i_tree;
        tree_count_subgroup_file[i_file]+=n_halos_tree_subgroup[i_tree];
      }
      n_halos_used+=tree_count_subgroup_file[i_file];
    }
  }
  while(i_tree<n_trees_subgroup-1){
    i_tree++;
    tree_hi_subgroup_file[n_files_subgroups-1]    =i_tree;
    tree_count_subgroup_file[n_files_subgroups-1]+=n_halos_tree_subgroup[i_tree];
  }

  // Determine how group trees will be distributed between files
  for(i_file=0,n_halos_used=0,i_tree=0;i_file<n_files_groups;i_file++,i_tree++){
    tree_lo_group_file[i_file]   =i_tree;
    tree_hi_group_file[i_file]   =i_tree;
    tree_count_group_file[i_file]=0;
    if(i_tree<n_trees_group){
      n_halos_target                =(n_halos_groups-n_halos_used)/(n_files_groups-i_file);
      tree_count_group_file[i_file]+=n_halos_tree_group[i_tree];
      while(tree_count_group_file[i_file]<n_halos_target && i_tree<n_trees_group-1){
        i_tree++;
        tree_hi_group_file[i_file]    =i_tree;
        tree_count_group_file[i_file]+=n_halos_tree_group[i_tree];
      }
      n_halos_used+=tree_count_group_file[i_file];
    }
  }
  while(i_tree<n_trees_group-1){
    i_tree++;
    tree_hi_group_file[n_files_groups-1]    =i_tree;
    tree_count_group_file[n_files_groups-1]+=n_halos_tree_group[i_tree];
  }

  tree_lo_subgroup_local   =tree_lo_subgroup_rank[SID.My_rank];
  tree_hi_subgroup_local   =tree_hi_subgroup_rank[SID.My_rank];
  n_trees_subgroup_local   =tree_hi_subgroup_local-tree_lo_subgroup_local+1;
  tree_count_subgroup_local=tree_count_subgroup_rank[SID.My_rank];
  tree_lo_group_local      =tree_lo_group_rank[SID.My_rank];
  tree_hi_group_local      =tree_hi_group_rank[SID.My_rank];
  n_trees_group_local      =tree_hi_group_local-tree_lo_group_local+1;
  tree_count_group_local   =tree_count_group_rank[SID.My_rank];

  // Report decomposition results
  if(n_files_subgroups>1){
    for(i_bin=0;i_bin<n_files_subgroups;i_bin++)
      SID_log("Subgroup file #%4d will store %5d subgroup trees (%8d halos in total)",SID_LOG_COMMENT,i_bin+1,tree_hi_subgroup_file[i_bin]-tree_lo_subgroup_file[i_bin]+1,tree_count_subgroup_file[i_bin]);    
  }
  if(n_files_groups>1){
    for(i_bin=0;i_bin<n_files_groups;i_bin++)
      SID_log("Group    file #%4d will store %5d group    trees (%8d halos in total)",SID_LOG_COMMENT,i_bin+1,tree_hi_group_file[i_bin]-tree_lo_group_file[i_bin]+1,tree_count_group_file[i_bin]);    
  }
  if(SID.n_proc>1){
    for(i_bin=0;i_bin<SID.n_proc;i_bin++)
      SID_log("Rank     #%4d will process %5d subgroup trees (%8d halos in total)",SID_LOG_COMMENT,i_bin+1,tree_hi_subgroup_rank[i_bin]-tree_lo_subgroup_rank[i_bin]+1,tree_count_subgroup_rank[i_bin]);
    for(i_bin=0;i_bin<SID.n_proc;i_bin++)
      SID_log("Rank     #%4d will process %5d group    trees (%8d halos in total)",SID_LOG_COMMENT,i_bin+1,tree_hi_group_rank[i_bin]-tree_lo_group_rank[i_bin]+1,tree_count_group_rank[i_bin]);
  }

  SID_log("(%d group and %d subgroup trees)...Done.",SID_LOG_CLOSE,n_trees_group,n_trees_subgroup);

  // VERTICAL TREE CONSTRUCTION STARTS HERE

  // Process subgroup trees (k_match==0) and then group trees (k_match==1)
  for(k_match=0;k_match<2;k_match++){
    switch(k_match){
      // Process subgroups
    case 0:
      n_halos_tree   =n_halos_tree_subgroup;
      n_trees_local  =n_trees_subgroup_local;
      tree_lo_rank   =tree_lo_subgroup_rank[SID.My_rank];
      tree_count_file=tree_count_subgroup_file;
      tree_lo_file   =tree_lo_subgroup_file;
      tree_hi_file   =tree_hi_subgroup_file;
      n_write        =n_files_subgroups;
      sprintf(group_text_prefix,"sub");
      break;
      // Process groups
    case 1:
      n_halos_tree   =n_halos_tree_group;
      n_trees_local  =n_trees_group_local;
      tree_lo_rank   =tree_lo_group_rank[SID.My_rank];
      tree_count_file=tree_count_group_file;
      tree_lo_file   =tree_lo_group_file;
      tree_hi_file   =tree_hi_group_file;
      n_write        =n_files_groups;
      sprintf(group_text_prefix,"");
      break;      
    }
    
    // Skip this iteration if the number of
    //   output files is set to <=0
    if(n_write<=0)
      continue;

    // Initialize tree arrays
    trees=(tree_info **)SID_malloc(sizeof(tree_info *)*n_trees_local);
    for(i_tree=0;i_tree<n_trees_local;i_tree++)
      init_tree(n_search,&trees[i_tree]);
    
    // Read matching files
    SID_log("Rendering %sgroup trees vertical...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
    if((*flag_clean))
      SID_log("(horizontal tree files will be removed)...",SID_LOG_CONTINUE);

    // Loop over all the horizontal tree files in order of decreasing snapshot number
    for(i_read=i_read_stop,i_file=i_read_stop,flag_init=TRUE;i_read>=i_read_start;i_read-=i_read_step,i_file--,flag_init=FALSE){
      halo_snap=(i_file)%n_search;

      // Only keep track of neighbour lists over the search range; initialize this snapshot's halo list here
      for(i_tree=0;i_tree<n_trees_local;i_tree++){
        trees[i_tree]->neighbour_halos[i_file%n_search]=NULL;
        trees[i_tree]->n_neighbours[i_file%n_search]   =0;
      }
      sprintf(filename_in,        "%s.trees_horizontal_%d",                    filename_root_out,i_read);
      sprintf(filename_properties,"%s.trees_horizontal_%sgroups_properties_%d",filename_root_out,group_text_prefix,i_read);

      // Open horizontal tree file
      fp_in=fopen(filename_in,"r");
      SID_fopen(filename_properties,"r",&fp_properties);

      // Read header
      grab_next_line(fp_in,&line,&line_length);
      grab_int(line,1,&n_groups);
      grab_int(line,2,&n_subgroups);
      
      // Read each group in turn
      for(i_group=0,i_subgroup=0,k_subgroup=0;i_group<n_groups;i_group++){
        j_read=1;
        grab_next_line(fp_in,&line,&line_length);
        grab_int(line,j_read++,&group_id);            // Unique id assigned to this group halo
        grab_int(line,j_read++,&group_descendant_id); // Unique id assigned to its descendant
        grab_int(line,j_read++,&group_tree_id);       // Unique id assigned to this group's tree
        grab_int(line,j_read++,&group_file_offset);   // Number of files back to the matched descendant
        grab_int(line,j_read++,&n_subgroups_group);   // Number of subgroups in this group
        group_tree_id=i_tree_group[group_tree_id];
        // If we are processing subgroup trees ...
        if(k_match==0){
          // Read each subgroup in turn
          for(i_subgroup=0,j_subgroup=4;i_subgroup<n_subgroups_group;i_subgroup++){
            grab_int(line,j_read++,&subgroup_id);              // Unique id assigned to this subgroup
            // Ignore negative ids
            if(subgroup_id>=0){
              grab_int(line,j_read++,&subgroup_descendant_id); // Unique id assigned to its descendant halo
              grab_int(line,j_read++,&subgroup_tree_id);       // Unique id assigned to this subgroup's tree
              grab_int(line,j_read++,&subgroup_file_offset);   // Number of files back to the matched descendant
              SID_fread_all(&properties,sizeof(halo_info),1,&fp_properties);
              subgroup_tree_id=i_tree_subgroup[subgroup_tree_id];
              i_tree          =subgroup_tree_id-tree_lo_rank;
              // If this subgroup belongs to a local tree ...
              if(i_tree>=0 && i_tree<n_trees_local){
                // ... create a new branch and add it to its tree ...
                if(subgroup_file_offset==0)
                  descendant_snap=-1;
                else
                  descendant_snap=(i_file+subgroup_file_offset)%n_search;
                add_node_to_tree(trees[i_tree],
                                 subgroup_id,
                                 group_id,
                                 subgroup_descendant_id,
                                 halo_snap,
                                 descendant_snap,
                                 &properties);
              }
            }
            // ... else if it is a dropped halo, skip 3 elements to the next subgroup
            else
              j_read+=3;
          }
        }
        // .. else, if we are processing group trees (again, ignore negative ids)
        else if(group_id>=0){
          SID_fread_all(&properties,sizeof(halo_info),1,&fp_properties);
          group_tree_id=i_tree_group[group_tree_id];
          i_tree       =group_tree_id-tree_lo_rank;
          // If this group belongs to a local tree ...
          if(i_tree>=0 && i_tree<n_trees_local){
            // ... if so, create a new branch and add it to the tree ...
            if(group_file_offset==0)
              descendant_snap=-1;
            else
              descendant_snap=(i_file+group_file_offset)%n_search;
            add_node_to_tree(trees[i_tree],
                             group_id,
                             group_id,
                             group_descendant_id,
                             halo_snap,
                             descendant_snap,                      
                             &properties);
          }
        }
      } // i_group
      fclose(fp_in);
      SID_fclose(&fp_properties);

      // Fix group halo ordering
      for(i_tree=0;i_tree<n_trees_local;i_tree++)
        assign_group_halo_order(trees[i_tree],halo_snap);      

      // If flag_clean=TRUE, then delete the input files used here.
      if((*flag_clean)==TRUE){
        if(k_match==1){
          sprintf(filename_in,"%s.trees_horizontal_%d",filename_root_out,i_read);
          remove(filename_in);
        }
        sprintf(filename_properties,"%s.trees_horizontal_%sgroups_properties_%d",filename_root_out,group_text_prefix,i_read);
        remove(filename_properties);
      }      
    } // i_read
    SID_log("Done.",SID_LOG_CLOSE);

    // Now that the trees are built, we can traverse them to assign halo ids.
    SID_log("Constructing %sgroup tree identities...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);

    // Make corrections to progenitor order ...
    SID_log("Correcting progenitor ordering...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
    for(i_tree=0;i_tree<n_trees_local;i_tree++){
      // Loop over each tree in the i_tree'th grouping.  They all have NULL descendants.
      current=trees[i_tree]->root;
      while(current!=NULL){
        if(current->descendant==NULL){
          progenitor_score=0;
          assign_progenitor_order_recursive(current,&progenitor_score,TREE_PROGENITOR_ORDER_DELUCIA);
        }
        current=current->next;
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);    

    //  ... then assign depth_first_indices ...
    SID_log("Assigning depth-first indices...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
    for(i_tree=0;i_tree<n_trees_local;i_tree++){
      depth_first_index=0;
      current=trees[i_tree]->root;
      while(current!=NULL){
        if(current->descendant==NULL)
          assign_depth_first_index_recursive(current,&depth_first_index);
        current=current->next;
      }
      if(n_halos_tree[i_tree]!=(int)depth_first_index)
        SID_trap_error("Halo count mismatch (ie. %d!=%d) in %sgroup tree %d",ERROR_LOGIC,n_halos_tree[i_tree],(int)depth_first_index,i_tree);
    }
    SID_log("Done.",SID_LOG_CLOSE);    

    //  ... and then assign unique ids
    SID_log("Assigning ids...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
    for(i_tree=0;i_tree<n_trees_local;i_tree++){
      current=trees[i_tree]->root;
      while(current!=NULL){
        if(current->descendant==NULL)
          assign_unique_ids_recursive(current,tree_lo_rank+i_tree);
        current=current->next;
      }
    }
    SID_log("Done.",SID_LOG_CLOSE);
     
    SID_log("Done.",SID_LOG_CLOSE);    
    
    SID_log("Writing %sgroup trees (structure size=%lld bytes)...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix,sizeof(halo_info));
    // Write the tree headers
    for(i_rank=0,i_write=0,i_tree=0,flag_write_init=TRUE;i_rank<SID.n_proc;i_rank++){
      if(SID.My_rank==i_rank){
        if(n_write==1){
          sprintf(filename_out,    "%s.%sgroup_trees",    filename_root_out,group_text_prefix);
          sprintf(filename_out_MBP,"%s.%sgroup_trees_MBP",filename_root_out,group_text_prefix);
        }
        else{
          sprintf(filename_out,    "%s.%sgroup_trees.%d",    filename_root_out,group_text_prefix,i_write);
          sprintf(filename_out_MBP,"%s.%sgroup_trees_MBP.%d",filename_root_out,group_text_prefix,i_write);
        }
        if(flag_write_init){
          SID_fopen(filename_out,    "w",&fp_out);
          SID_fopen(filename_out_MBP,"w",&fp_out_MBP);
          n_trees_file=tree_hi_file[i_write]-tree_lo_file[i_write]+1;
          SID_fwrite(&n_trees_file,sizeof(int),1,&fp_out);
          SID_fwrite(&n_trees_file,sizeof(int),1,&fp_out_MBP);
          SID_fwrite(&(tree_count_file[i_write]),sizeof(int),1,&fp_out);
          SID_fwrite(&(tree_count_file[i_write]),sizeof(int),1,&fp_out_MBP);
        }
        else{
          SID_fopen(filename_out,    "a",&fp_out);
          SID_fopen(filename_out_MBP,"a",&fp_out_MBP);
        }
        // Write all of the trees on this rank
        for(j_tree=0;j_tree<n_trees_local;i_tree++,j_tree++){
            
          // Write halo counts for trees
          SID_fwrite(&(n_halos_tree[i_tree]),sizeof(int),1,&fp_out);
          SID_fwrite(&(n_halos_tree[i_tree]),sizeof(int),1,&fp_out_MBP);
          flag_write_init=FALSE;

          // Move to the next file if this one is done
          if(i_tree==tree_hi_file[i_write]){
            i_write++;
            if(i_write<n_write){
              flag_write_init=TRUE;
              SID_fclose(&fp_out);
              SID_fclose(&fp_out_MBP);
              sprintf(filename_out,    "%s.%sgroup_trees.%d",    filename_root_out,group_text_prefix,i_write);
              sprintf(filename_out_MBP,"%s.%sgroup_trees_MBP.%d",filename_root_out,group_text_prefix,i_write);
              SID_fopen(filename_out,    "w",&fp_out);
              SID_fopen(filename_out_MBP,"w",&fp_out_MBP);
              n_trees_file=tree_hi_file[i_write]-tree_lo_file[i_write]+1;
              SID_fwrite(&n_trees_file,sizeof(int),1,&fp_out);
              SID_fwrite(&n_trees_file,sizeof(int),1,&fp_out_MBP);
              SID_fwrite(&(tree_count_file[i_write]),sizeof(int),1,&fp_out);
              SID_fwrite(&(tree_count_file[i_write]),sizeof(int),1,&fp_out_MBP);
            }
          }
        }
        SID_fclose(&fp_out);
        SID_fclose(&fp_out_MBP);
      }
      // Update the other ranks on i_rank's progress
#ifdef USE_MPI
      MPI_Bcast(&i_write,        1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
      MPI_Bcast(&i_tree,         1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
      MPI_Bcast(&flag_write_init,1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
#endif
    }
    
    // Write the tree halos
    for(i_rank=0,i_write=0,i_tree=0,flag_write_init=TRUE;i_rank<SID.n_proc;i_rank++){
      if(SID.My_rank==i_rank){
        if(n_write==1){
          sprintf(filename_out,    "%s.%sgroup_trees",    filename_root_out,group_text_prefix);
          sprintf(filename_out_MBP,"%s.%sgroup_trees_MBP",filename_root_out,group_text_prefix);
        }
        else{
          sprintf(filename_out,    "%s.%sgroup_trees.%d",    filename_root_out,group_text_prefix,i_write);
          sprintf(filename_out_MBP,"%s.%sgroup_trees_MBP.%d",filename_root_out,group_text_prefix,i_write);
        }
        SID_fopen(filename_out,    "a",&fp_out);
        SID_fopen(filename_out_MBP,"a",&fp_out_MBP);
        // Write all of the trees on this rank
        for(j_tree=0,n_halos_written=0;j_tree<n_trees_local;i_tree++,j_tree++){
          // Loop over each tree in this i_tree'th grouping.  They all have NULL descendants.
          current=trees[j_tree]->root;
          while(current!=NULL){
            if(current->descendant==NULL)
              n_halos_written+=write_tree_vertical_halos_recursive(current,&fp_out,&fp_out_MBP);
            current=current->next;
          }
          flag_write_init=FALSE;
          // Move to the next file if this one is done
          if(i_tree==tree_hi_file[i_write]){
            i_write++;
            if(i_write<n_write){
              flag_write_init=TRUE;
              SID_fclose(&fp_out);
              SID_fclose(&fp_out_MBP);
              sprintf(filename_out,    "%s.%sgroup_trees.%d",    filename_root_out,group_text_prefix,i_write);
              sprintf(filename_out_MBP,"%s.%sgroup_trees_MBP.%d",filename_root_out,group_text_prefix,i_write);
              SID_fopen(filename_out,    "a",&fp_out);
              SID_fopen(filename_out_MBP,"a",&fp_out_MBP);
            }
          }
        }
        SID_fclose(&fp_out);
        SID_fclose(&fp_out_MBP);
      }
      // Update the other ranks on i_rank's progress
#ifdef USE_MPI
      MPI_Bcast(&i_write,        1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
      MPI_Bcast(&i_tree,         1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
      MPI_Bcast(&flag_write_init,1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
      MPI_Bcast(&n_halos_written,1,MPI_INTEGER,i_rank,MPI_COMM_WORLD);
#endif
    }
    SID_log("(%d trees and %d halos written)...Done.",SID_LOG_CLOSE,i_tree,n_halos_written);
    // Free trees
    for(i_tree=0;i_tree<n_trees_local;i_tree++)
      free_tree(&trees[i_tree]);
    SID_free((void **)&trees);
  } // k_match
  
  // Clean-up
  SID_free((void **)&line);
  SID_free((void **)&i_tree_group);
  SID_free((void **)&i_tree_subgroup);
  SID_free((void **)&n_halos_tree_group);
  SID_free((void **)&n_halos_tree_subgroup);
  SID_free((void **)&tree_lo_group_rank);
  SID_free((void **)&tree_hi_group_rank);
  SID_free((void **)&tree_count_group_rank);
  SID_free((void **)&tree_lo_subgroup_rank);
  SID_free((void **)&tree_hi_subgroup_rank);
  SID_free((void **)&tree_count_subgroup_rank);
  SID_free((void **)&tree_lo_group_file);
  SID_free((void **)&tree_hi_group_file);
  SID_free((void **)&tree_count_group_file);
  SID_free((void **)&tree_lo_subgroup_file);
  SID_free((void **)&tree_hi_subgroup_file);
  SID_free((void **)&tree_count_subgroup_file);

  SID_log("Done.",SID_LOG_CLOSE);
}
