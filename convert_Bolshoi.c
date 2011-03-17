#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

#define N_UPIDS_MAX 100000

void add_uberparent_to_list(int tree_id,int *tree_ids,int *n_tree_ids);
void add_uberparent_to_list(int tree_id,int *tree_ids,int *n_tree_ids){
  int flag_found;
  int i_tree_id;
  for(i_tree_id=0,flag_found=FALSE;i_tree_id<(*n_tree_ids) && !flag_found;i_tree_id++){
    if(tree_id==tree_ids[i_tree_id]){
      flag_found=TRUE;
      break;
    }
  }
  if(!flag_found)
    tree_ids[(*n_tree_ids)++]=tree_id;
}

void add_scale_to_list(float scale,float *scales,int *n_scales);
void add_scale_to_list(float scale,float *scales,int *n_scales){
  int flag_found;
  int i_scale;
  for(i_scale=0,flag_found=FALSE;i_scale<(*n_scales) && !flag_found;i_scale++){
    if(scale==scales[i_scale]){
      flag_found=TRUE;
      break;
    }
  }
  if(!flag_found)
    scales[(*n_scales)++]=scale;
}

int scale_to_snap(float scale,float *scales,int n_scales);
int scale_to_snap(float scale,float *scales,int n_scales){
  int flag_found;
  int i_scale;
  for(i_scale=0,flag_found=FALSE;i_scale<n_scales && !flag_found;i_scale++){
    if(scale==scales[i_scale]){
      flag_found=TRUE;
      break;
    }
  }
  if(!flag_found){
    SID_log("List of available scale factors:",SID_LOG_OPEN);
    for(i_scale=0;i_scale<n_scales;i_scale++)
      SID_log("%f",SID_LOG_COMMENT,scales[i_scale]);
    SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
    SID_trap_error("Could not find scale=%f in list of scales!",ERROR_LOGIC,scale);
  }
  return(i_scale);
}

float propagate_spins_recursive(tree_node_info *tree,RNG_info *RNG);
float propagate_spins_recursive(tree_node_info *tree,RNG_info *RNG){
  tree_node_info *current;
  float           spin_propagate;
  float           spin;
  float           V_vir;

  // This is what gets returned if there are no progenitors
  spin_propagate=tree->halo.spin[0];
  // Walk the tree
  current=tree->progenitor_first;
  while(current!=NULL){
    spin=propagate_spins_recursive(current,RNG);
    // If this is the first progenitor, then it's spin gets propagated up the tree.
    if(current==tree->progenitor_first)
      spin_propagate=spin;
    current=current->progenitor_next;
  }

  // Convert lambda -> j_specific here.
  V_vir             =1e-3*sqrt((double)(G_NEWTON*(double)(tree->halo.M_vir)*M_SOL*1e10/((double)(tree->halo.R_vir)*M_PER_MPC)));
  tree->halo.spin[0]=spin_propagate*sqrt(2.)*tree->halo.R_vir*V_vir;
  return(spin_propagate);
}

int main(int argc, char *argv[]){
  char        filename_in_root[256];
  char        filename_out_root[256];
  char        filename_in[256];
  char        filename_out[256];
  char        filename_snaps[256];
  char       *line=NULL;
  size_t      line_length=0;
  FILE       *fp_in;
  FILE       *fp_out;
  FILE       *fp_snaps;
  char        comment_char[2];
  char        first_char[2];
  int   n_halos;
  int   n_trees_out,n_trees_in;
  int  *n_halos_tree;
  int  *n_trees_tree;
  int   i_halo,j_halo;
  int   i_tree;
  int   j_tree;
  int   n_tree_ids;
  int   tree_ids[N_UPIDS_MAX];
  float scales[250];
  float scales_check[250];
  int   n_halos_scale[250];
  int   n_scales,i_scale,n_scales_check;
  int   flag_found;
  float descendant_scale;
  float halo_scale;
  int   halo_id;
  int   group_id;
  int   search_id;
  int   flag_init=TRUE;
  int   n_halos_written;
  int   descendant_id;
  int   num_prog;
  int   parent_id;
  int   uberparent_id;
  int   descendant_parent_id;
  int   merger_type;
  float M_vir;
  float M_trunc;
  float R_vir;
  float V_vir;
  float r_s;
  float R_trunc;
  int   flag_MMP;
  float LMM_scale;
  float V_max;
  float x;
  float y;
  float z;
  float v_x;
  float v_y;
  float v_z;
  float sigma_v;
  float spin_x;
  float spin_y;
  float spin_z;
  float m_p;
  int   n_particles;
  int   n_header_lines;
  int   i_header;

  int   i_write_start,i_write_stop;
  int   i_report,i_tree_next_report;

  int   halo_snap;
  int   descendant_snap;
  int   LMM_snap;
  int   progenitor_score;

  tree_info      **trees;
  halo_info        properties;
  tree_node_info  *current;
  tree_node_info  *next;

  int               tree_mode;
  int               progenitor_mode;
  int               depth_first_index;
  int               flag_clean;
  int               i_x,i_y,i_z,i_write,n_write;
  char              filename_root_in[256];
  char              filename_root_out[256];
  RNG_info          RNG;
  int               seed=102873;
  
  SID_init(&argc,&argv,NULL);

  SID_log("Converting Bolshoi trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Fetch user inputs
  m_p              =1.35e8;
  n_write          =5*5*5;
  n_header_lines   =8;

  i_write_start  =atoi(argv[1]);
  i_write_stop   =atoi(argv[2]);

  //progenitor_mode=TREE_PROGENITOR_ORDER_DELUCIA;
  progenitor_mode=TREE_PROGENITOR_ORDER_DEFAULT;
  tree_mode      =0; // join-on

  init_RNG(&seed,&RNG,RNG_DEFAULT);

  // Loop over all raw files
  for(i_x=0,i_write=0;i_x<5 && i_write<n_write;i_x++){
    for(i_y=0;i_y<5 && i_write<n_write;i_y++){
      for(i_z=0;i_z<5 && i_write<n_write;i_z++,i_write++){
        if(i_write>=i_write_start && i_write<=i_write_stop){
          SID_log("Processing file %d out of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_write+1,n_write);

          // Open file
          sprintf(filename_in,"/nfs/dset/shrek071/millenium/bolshoi/raw/tree_%d_%d_%d.dat",i_x,i_y,i_z);
          SID_log("Open file {%s}...",SID_LOG_OPEN,filename_in);
          fp_in=fopen(filename_in,"r");
          SID_log("Done.",SID_LOG_CLOSE);

          // Skip comment lines at the top
          SID_log("Skip header comments...",SID_LOG_OPEN);
          for(i_header=0;i_header<n_header_lines;i_header++)
            grab_next_line(fp_in,&line,&line_length);
          SID_log("Done.",SID_LOG_CLOSE);

          // Read number of iso-trees
          SID_log("Read number of iso-trees...",SID_LOG_OPEN);
          grab_next_line(fp_in,&line,&line_length);
          grab_int(line,1,&n_trees_in);   
          SID_log("{%d}...Done.",SID_LOG_CLOSE,n_trees_in);

          // Read list of expansion factors
          SID_log("Reading scale list...",SID_LOG_OPEN);
          sprintf(filename_snaps,"/nfs/dset/shrek071/millenium/bolshoi/treedata/Bolshoi.a_list",i_write);
          fp_snaps=fopen(filename_snaps,"r");
          n_scales=count_lines_data(fp_snaps);
          for(i_scale=0;i_scale<n_scales;i_scale++){
            grab_next_line_data(fp_snaps,&line,&line_length);
            grab_float(line,1,&(scales[i_scale]));
          }
          fclose(fp_snaps);
          SID_log("(%d found)...",SID_LOG_CONTINUE,n_scales);
          SID_log("Done.",SID_LOG_CLOSE);

          // Ensure that the scales in all files are the same
          if(i_write==i_write_start){
            n_scales_check=n_scales;
            for(i_scale=0;i_scale<n_scales;i_scale++)
              scales_check[i_scale]=scales[i_scale];
          }
          else{
            if(n_scales!=n_scales_check)
              SID_trap_error("n_scales mismatch (ie. %d!=%d)",ERROR_LOGIC,n_scales,n_scales_check);
            for(i_scale=0;i_scale<n_scales;i_scale++){
              if(scales[i_scale]!=scales_check[i_scale])
                SID_trap_error("scale mismatch (ie. %lf!=%lf)",ERROR_LOGIC,scales[i_scale],scales_check[i_scale]);
            }          
          }        

          // Count input-trees and halos and determine output-tree ownership
          SID_log("Counting trees and halos...",SID_LOG_OPEN|SID_LOG_TIMER);
          grab_next_line(fp_in,&line,&line_length); // First iso-tree header line
          i_tree  =0;
          n_halos =0;
          n_tree_ids =0;
          // Read each tree in turn
          while(!feof(fp_in)){
            // Process iso-tree halos
            i_halo=0;
            do{
              grab_next_line(fp_in,&line,&line_length);
              if(strlen(line)!=0){ // Ignore blank lines
                // Increment counters (if this isn't a new iso-tree header line)
                if(!check_comment(line)){
                    if(i_halo==0){
                      grab_int(line,7,&uberparent_id);
                      grab_int(line,6,&parent_id);
                      grab_int(line,2,&halo_id);
                      if(uberparent_id>=0)
                        group_id=uberparent_id;
                      else if(parent_id>=0)
                        group_id=parent_id;
                      else
                        group_id=halo_id;
                      if(tree_mode==1)
                        group_id=halo_id;
                      else if(tree_mode==2)
                        group_id=0;
                      add_uberparent_to_list(group_id,tree_ids,&n_tree_ids);
                      if(n_tree_ids>N_UPIDS_MAX)
                        SID_trap_error("Increase N_UPIDS_MAX!",ERROR_LOGIC);
                    }
                    n_halos++;
                    i_halo++;
                }
              }
            } while(!check_comment(line));
            i_tree++;
          }
          if(i_tree!=n_trees_in)
            SID_trap_error("The number of trees read (%d) does not equal the number listed in the header (%d)!",ERROR_LOGIC,i_tree,n_trees_in);
          merge_sort(tree_ids,(size_t)n_tree_ids,NULL,SID_INT,SORT_INPLACE_ONLY,TRUE);
          SID_log("(%d trees and %d halos found; %d output-trees will result)...Done.",SID_LOG_CLOSE,n_trees_in,n_halos,n_tree_ids);

          // Allocate trees
          SID_log("Initializing trees...",SID_LOG_OPEN|SID_LOG_TIMER);
          n_trees_out =n_tree_ids;
          trees       =(tree_info **)SID_malloc(sizeof(tree_info *)*n_trees_out);
          n_halos_tree=(int        *)SID_malloc(sizeof(int)*n_trees_out);
          n_trees_tree=(int        *)SID_malloc(sizeof(int)*n_trees_out);
          for(i_tree=0;i_tree<n_trees_out;i_tree++){
            init_tree(n_scales,&(trees[i_tree]));
            n_halos_tree[i_tree]=0;
            n_trees_tree[i_tree]=0;
          }
          SID_log("Done.",SID_LOG_CLOSE);

          // Rewind and skip header
          rewind(fp_in);
          for(i_header=0;i_header<n_header_lines;i_header++)
            grab_next_line(fp_in,&line,&line_length);
          grab_next_line(fp_in,&line,&line_length);

          // FILE CONVERSION STARTS HERE
          SID_log("Processing trees...",SID_LOG_OPEN|SID_LOG_TIMER);
          grab_next_line(fp_in,&line,&line_length);
          i_tree  =0;
          i_halo  =0;
          i_report=0;
          i_tree_next_report=n_trees_in/10;
          while(!feof(fp_in)){
            j_halo=0;
            do{
              grab_next_line(fp_in,&line,&line_length);
              if(strlen(line)!=0 && !check_comment(line)){ // Ignore blank lines

                // Parse this halo's line
                /*
                  parse_line(line,17,
                  1, &halo_scale,          SID_FLOAT,
                  2, &halo_id,             SID_INT,
                  3, &descendant_scale,    SID_FLOAT,
                  4, &descendant_id,       SID_INT,
                  5, &num_prog,            SID_INT,
                  6, &parent_id,           SID_INT,
                  7, &uberparent_id,       SID_INT,
                  8, &descendant_parent_id,SID_INT,
                  9, &merger_type,         SID_INT,
                  10,&M_vir,               SID_FLOAT,
                  11,&M_trunc,             SID_FLOAT,
                  12,&R_vir,               SID_FLOAT,
                  13,&r_s,                 SID_FLOAT,
                  14,&R_trunc,             SID_FLOAT,
                  15,&flag_MMP,            SID_FLOAT,
                  16,&LMM_scale,           SID_FLOAT,
                  17,&V_max,               SID_FLOAT);
                */
                grab_float(line,1, &halo_scale);
                grab_int(  line,2, &halo_id);
                grab_float(line,3, &descendant_scale);
                grab_int(  line,4, &descendant_id);
                grab_int(  line,6, &parent_id);
                grab_int(  line,7, &uberparent_id);
                grab_float(line,11,&M_vir);
                grab_float(line,12,&R_vir);
                grab_float(line,14,&sigma_v);
                grab_float(line,17,&V_max);
                grab_float(line,18,&x);
                grab_float(line,19,&y);
                grab_float(line,20,&z);
                grab_float(line,21,&v_x);
                grab_float(line,22,&v_y);
                grab_float(line,23,&v_z);

                // Determine which output-tree this input-tree belongs to
                if(j_halo==0){
                  if(uberparent_id>=0)
                    search_id=uberparent_id;
                  else if(parent_id>=0)
                    search_id=parent_id;
                  else
                    search_id=halo_id;
                  if(tree_mode==1)
                    j_tree=i_tree;
                  else if(tree_mode==2)
                    j_tree=0;
                  else{
                    j_tree=find_index_int(tree_ids,search_id,n_tree_ids,NULL);
                    if(tree_ids[j_tree]!=search_id)
                      SID_trap_error("Arg! Couldn't find tree_id=%d in the list",ERROR_LOGIC,search_id);
                  }
                  n_trees_tree[j_tree]++;
                }

                // Set the group id
                if(uberparent_id>=0)
                  group_id=uberparent_id;
                else if(parent_id>=0)
                  group_id=parent_id;
                else
                  group_id=halo_id;

                // Generate lambdas
                V_vir  =1e-3*sqrt(G_NEWTON*M_vir*M_SOL/(R_vir*M_PER_KPC));
                spin_x =(float)random_lognormal(&RNG,0.045,0.5);              // Constants taken from Bullock et al '01 (also see Vitvitska '02)
                spin_y =0.;
                spin_z =0.;

                // Convert scales to snapshot indices
                if(descendant_scale==0.)
                  descendant_snap=-1;
                else
                  descendant_snap=scale_to_snap(descendant_scale,scales,n_scales);
                halo_snap=scale_to_snap(halo_scale,scales,n_scales);

                // Set the needed halo properties
                properties.n_particles  =(int)(M_vir/m_p);
                properties.M_vir        =M_vir/1e10;
                properties.R_vir        =R_vir/1e3;
                properties.pos[0]       =x;
                properties.pos[1]       =y;
                properties.pos[2]       =z;
                properties.vel[0]       =v_x;
                properties.vel[1]       =v_y;
                properties.vel[2]       =v_z;
                properties.sigma_v      =sigma_v;
                properties.v_max        =V_max;
                properties.spin[0]      =spin_x;
                properties.spin[1]      =spin_y;
                properties.spin[2]      =spin_z;
                properties.most_bound_id=0;
                properties.snap_num     =halo_snap;
                properties.halo_index   =i_halo;
                properties.halo_id      =halo_id;
                properties.group_id     =group_id;

                // Add this halo to its tree
                add_node_to_tree(trees[j_tree],
                                 halo_id,
                                 group_id,
                                 descendant_id,
                                 halo_snap,
                                 descendant_snap,
                                 &properties);
                i_halo++;
                j_halo++;
                n_halos_tree[j_tree]++;
              }
            } while(!check_comment(line));

            // Write a status message
            if(i_tree==i_tree_next_report){
              i_report++;
              SID_log("%3d%% complete.",SID_LOG_COMMENT|SID_LOG_TIMER,10*i_report);
              i_tree_next_report=MIN(n_trees_in,n_trees_in*(i_report+1)/10);
            } 
            i_tree++;
          }
          fclose(fp_in);
          SID_log("Done.",SID_LOG_CLOSE);

          // Write trees
          sprintf(filename_out,"/nfs/dset/shrek071/millenium/bolshoi/treedata_wip2/trees_%d.%d",n_scales-1,i_write);
          SID_log("Writing output file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_out);
          fp_out =fopen(filename_out,"w");

          // ... write header ...
          SID_log("Writing header...",SID_LOG_OPEN|SID_LOG_TIMER);
          fwrite(&n_trees_out,sizeof(int),1,fp_out);
          fwrite(&n_halos,    sizeof(int),1,fp_out);
          for(i_tree=0;i_tree<n_trees_out;i_tree++)
            fwrite(&(n_halos_tree[i_tree]),sizeof(int),1,fp_out);
          SID_log("Done.",SID_LOG_CLOSE);

          // ... construct IDs ...
          SID_log("Constructing IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
          for(i_tree=0;i_tree<n_trees_out;i_tree++){
            // ... correct group halo ordering ...
            for(i_scale=0;i_scale<n_scales;i_scale++)
              assign_group_halo_order(trees[i_tree],i_scale);      
            // ... correct progenitor ordering ...
            current=trees[i_tree]->root;
            while(current!=NULL){
              if(current->descendant==NULL){
                progenitor_score=0;
                assign_progenitor_order_recursive(current,&progenitor_score,progenitor_mode);
              }
              current=current->next;
            }
            // ... assign depth-first-indices ...
            depth_first_index=0;
            current=trees[i_tree]->root;
            while(current!=NULL){
              if(current->descendant==NULL)
                assign_depth_first_index_recursive(current,&depth_first_index);
              current=current->next;
            }
            if(depth_first_index!=n_halos_tree[i_tree])
              SID_trap_error("DFI != n_halos (i.e. %d!=%d)",ERROR_LOGIC,depth_first_index,n_halos_tree[i_tree]);
            // ... assign ids ...
            current=trees[i_tree]->root;
            while(current!=NULL){
              if(current->descendant==NULL)
                assign_unique_ids_recursive(current,i_tree);
              current=current->next;
            }
          }
          SID_log("Done.",SID_LOG_CLOSE);

          // ... construct spins (needs to be done after IDs are assigned)...
          SID_log("Propagating spins...",SID_LOG_OPEN|SID_LOG_TIMER);
          for(i_tree=0;i_tree<n_trees_out;i_tree++){
            current=trees[i_tree]->root;
            while(current!=NULL){
              if(current->descendant==NULL)
                propagate_spins_recursive(current,&RNG);
              current=current->next;
            }
          }
          SID_log("Done.",SID_LOG_CLOSE);

          // ... write trees ...
          SID_log("Writing trees...",SID_LOG_OPEN|SID_LOG_TIMER);
          for(i_tree=0;i_tree<n_trees_out;i_tree++){
            n_halos_written=0;
            current=trees[i_tree]->root;
            while(current!=NULL){
              if(current->descendant==NULL)
                n_halos_written+=write_tree_vertical_halos_recursive(current,&fp_out,NULL);
              current=current->next;
            }
            if(n_halos_written!=n_halos_tree[i_tree])
              SID_trap_error("Number of halos written is not right (i.e. %d!=%d) for tree #%d",ERROR_LOGIC,n_halos_written,n_halos_tree[i_tree],i_tree);
          }
          SID_log("Done.",SID_LOG_CLOSE);
          SID_log("Done.",SID_LOG_CLOSE);

          // Free trees
          SID_log("Freeing trees...",SID_LOG_OPEN|SID_LOG_TIMER);
          for(i_tree=0;i_tree<n_trees_out;i_tree++)
            free_tree(&(trees[i_tree]));
          SID_log("Done.",SID_LOG_CLOSE);

          // Print some stats for this input file
          merge_sort(n_trees_tree,(size_t)n_trees_out,NULL,SID_INT,SORT_INPLACE_ONLY,TRUE);
          SID_log("Min # of input_trees to a tree: %d",SID_LOG_COMMENT,n_trees_tree[0]);
          SID_log("Med # of input_trees to a tree: %d",SID_LOG_COMMENT,n_trees_tree[n_trees_out/2]);
          SID_log("Max # of input_trees to a tree: %d",SID_LOG_COMMENT,n_trees_tree[n_trees_out-1]);

          // Clean-up
          SID_free(SID_FARG trees);
          SID_free(SID_FARG n_halos_tree);
          SID_free(SID_FARG n_trees_tree);
          fclose(fp_out);

          SID_log("Done.",SID_LOG_CLOSE);
        }
      }
    }
  }
  free_RNG(&RNG);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
