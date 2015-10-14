#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

#define _FILE_OFFSET_BITS 64

int main(int argc, char *argv[]){
  char        filename_root_in[MAX_FILENAME_LENGTH];
  char        filename_in[MAX_FILENAME_LENGTH];
  char        filename_out[MAX_FILENAME_LENGTH];
  char        group_text_prefix[8];
  int         grid_size;
  double      box_size;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  strcpy(filename_root_in, argv[1]);
  grid_size=atoi(argv[2]);
  box_size =(double)atof(argv[3]);

  // Determine the filename root of the in/out directories
  char filename_name_in[MAX_FILENAME_LENGTH];
  strcpy(filename_name_in,filename_root_in);
  strip_path(filename_name_in);

  SID_log("Re-writing trees {%s} to a %dX%dX%d grid...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in,grid_size,grid_size,grid_size);

  // Process subgroups and then groups
  int k_match;
  for(k_match=0;k_match<2;k_match++){
     switch(k_match){
        case 0:
           sprintf(group_text_prefix,"sub");
           break;
        case 1:
           sprintf(group_text_prefix,"");
           break;
     }
     SID_log("Processing %sgroup trees...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);

     // Counting input files
     int   n_in=0;
     int   flag_exists=TRUE;
     FILE *fp_test;
     int   n_trees_total=0;
     while(flag_exists){
        sprintf(filename_in,"%s/vertical/%sgroup_trees_%03d.dat",filename_root_in,group_text_prefix,n_in);
        fp_test=fopen(filename_in,"r");
        if(fp_test==NULL)
           flag_exists=FALSE;
        else{
           int n_trees_file;
           fread(&n_trees_file,sizeof(int),1,fp_test);
           n_trees_total+=n_trees_file;
           fclose(fp_test);
           n_in++;
        }
     }

     if(n_in==0)
        SID_trap_error("No %sgroup files found.",ERROR_LOGIC,group_text_prefix);
     SID_log("%d input files found.",SID_LOG_COMMENT,n_in);

     // Create some arrays of counters
     int   n_files_out;
     int  *n_halos_file_out;
     int  *n_trees_file_out;
     n_files_out     =grid_size*grid_size*grid_size;
     n_halos_file_out=(int  *)SID_calloc(sizeof(int)*n_files_out);
     n_trees_file_out=(int  *)SID_calloc(sizeof(int)*n_files_out);

     // Count the number of input trees and halos that will end-up in each output file
     SID_log("Counting the trees and halos that will be in each output file...",SID_LOG_OPEN|SID_LOG_TIMER,n_in);
     int        i_file;
     int        k_tree;
     FILE      *fp_in;
     int        n_halos_biggest_tree;
     int        i_report,next_report;
     halo_properties_SAGE_info *halo;
     halo=(halo_properties_SAGE_info *)SID_calloc(sizeof(halo_properties_SAGE_info));
     n_halos_biggest_tree=0;
     i_report            =1;
     next_report         =(int)((float)i_report*0.1*(float)n_trees_total);
     for(i_file=0,k_tree=0;i_file<n_in;i_file++){
        // Open file and read header
        int  n_trees_file;
        int  n_halos_file;
        int  i_tree;
        int *n_halos_tree;
        sprintf(filename_in, "%s/vertical/%sgroup_trees_%03d.dat",filename_root_in,group_text_prefix,i_file);
        fp_in=fopen(filename_in,"r");
        // Read the input file's header
        fread(&n_trees_file,sizeof(int),1,fp_in);
        fread(&n_halos_file,sizeof(int),1,fp_in);
        n_halos_tree=(int *)SID_malloc(sizeof(int)*n_trees_file);
        fread(n_halos_tree,sizeof(int),n_trees_file,fp_in);
        // Read the first halo of each tree to decide where it goes
        for(i_tree=0;i_tree<n_trees_file;i_tree++,k_tree++){
           int i_x,i_y,i_z,i_g;
           // Read the first halo and skip the rest
           fread(halo,sizeof(halo_properties_SAGE_info),1,fp_in);
           fseeko(fp_in,(size_t)(n_halos_tree[i_tree]-1)*sizeof(halo_properties_SAGE_info),SEEK_CUR);
           // Decide which grid cell it belongs to
           i_x=(int)((float)grid_size*(halo->pos[0]/box_size));i_x=MIN(i_x,grid_size-1);
           i_y=(int)((float)grid_size*(halo->pos[1]/box_size));i_y=MIN(i_y,grid_size-1);
           i_z=(int)((float)grid_size*(halo->pos[2]/box_size));i_z=MIN(i_z,grid_size-1);
           i_g=(i_z*grid_size+i_y)*grid_size+i_x;
           // Calculate some statistics
           n_halos_file_out[i_g]+=n_halos_tree[i_tree];
           n_trees_file_out[i_g]++;
           if(n_halos_biggest_tree<n_halos_tree[i_tree]) 
              n_halos_biggest_tree=n_halos_tree[i_tree];
           // Status message
           if(k_tree==next_report){
              SID_log("%3d%% complete.",SID_LOG_COMMENT|SID_LOG_TIMER,i_report*10);
              i_report++;
              if(i_report==10)
                 next_report=n_trees_total-1;
              else
                 next_report=(int)((float)i_report*0.1*(float)n_trees_total);
           }
        }
        fclose(fp_in);
        SID_free(SID_FARG n_halos_tree);
     }
     SID_log("Done.",SID_LOG_CLOSE);

     // Repeat to build-up a list of tree sizes needed for the header of each output file
     SID_log("Determining lists of tree sizes for each output file...",SID_LOG_OPEN|SID_LOG_TIMER,n_in);
     // Create the arrays containing the tree sizes which will be written to the output file's headers
     int  **tree_size;
     tree_size=(int **)SID_malloc(sizeof(int *)*n_files_out);
     for(i_file=0;i_file<n_files_out;i_file++){
        tree_size[i_file]=(int *)SID_calloc(sizeof(int)*n_trees_file_out[i_file]);
        n_trees_file_out[i_file]=0; // re-set to zero so we can use it as a counter
     }
     i_report   =1;
     next_report=(int)((float)i_report*0.1*(float)n_trees_total);
     for(i_file=0,k_tree=0;i_file<n_in;i_file++){
        // Open input file
        int  n_trees_file;
        int  n_halos_file;
        int  i_tree;
        int *n_halos_tree;
        sprintf(filename_in, "%s/vertical/%sgroup_trees_%03d.dat",filename_root_in,group_text_prefix,i_file);
        fp_in=fopen(filename_in,"r");
        // Read the input file's header
        fread(&n_trees_file,sizeof(int),1,fp_in);
        fread(&n_halos_file,sizeof(int),1,fp_in);
        n_halos_tree=(int *)SID_malloc(sizeof(int)*n_trees_file);
        fread(n_halos_tree,sizeof(int),n_trees_file,fp_in);
        for(i_tree=0;i_tree<n_trees_file;i_tree++,k_tree++){
           if(n_halos_tree[i_tree]>0){
              int i_x,i_y,i_z,i_g;
              // Read the first halo and skip the rest
              fread(halo,sizeof(halo_properties_SAGE_info),1,fp_in);
              fseeko(fp_in,(size_t)(n_halos_tree[i_tree]-1)*sizeof(halo_properties_SAGE_info),SEEK_CUR);
              // Decide which grid cell it belongs to
              i_x=(int)((float)grid_size*(halo->pos[0]/box_size));i_x=MIN(i_x,grid_size-1);
              i_y=(int)((float)grid_size*(halo->pos[1]/box_size));i_y=MIN(i_y,grid_size-1);
              i_z=(int)((float)grid_size*(halo->pos[2]/box_size));i_z=MIN(i_z,grid_size-1);
              i_g=(i_z*grid_size+i_y)*grid_size+i_x;
              // Add this tree to the right grid cell's tree size list
              tree_size[i_g][n_trees_file_out[i_g]++]=n_halos_tree[i_tree];
              // Status message
              if(k_tree==next_report){
                 SID_log("%3d%% complete.",SID_LOG_COMMENT|SID_LOG_TIMER,i_report*10);
                 i_report++;
                 if(i_report==10)
                    next_report=n_trees_total-1;
                 else
                    next_report=(int)((float)i_report*0.1*(float)n_trees_total);
              }
           }
        }
        fclose(fp_in);
        SID_free(SID_FARG n_halos_tree);
     }
     SID_free(SID_FARG halo);
     SID_log("Done.",SID_LOG_CLOSE);

     // Write the headers for each output file
     FILE  *fp_out;
     sprintf(filename_out,"%s/vertical_grid/",
             filename_root_in);
     SID_log("Creating output directory {%s}...",SID_LOG_OPEN,filename_out);
     mkdir(filename_out,02755);
     SID_log("Done.",SID_LOG_CLOSE);
     SID_log("Writing output file headers...",SID_LOG_OPEN,filename_out);
     for(i_file=0;i_file<n_files_out;i_file++){
        sprintf(filename_out,"%s/vertical_grid/%sgroup_trees_%03d.dat",filename_root_in,group_text_prefix,i_file);
        fp_out=fopen(filename_out,"w");
        fwrite(&(n_trees_file_out[i_file]),sizeof(int),1,                       fp_out);
        fwrite(&(n_halos_file_out[i_file]),sizeof(int),1,                       fp_out);
        fwrite(tree_size[i_file],          sizeof(int),n_trees_file_out[i_file],fp_out);
        fclose(fp_out);
        n_trees_file_out[i_file]=0; // re-set to zero so we can use it as a counter
     }
     SID_log("Done.",SID_LOG_CLOSE);

     // Read each tree and write it to the appropriate file
     SID_log("Distributing input trees to the grid of output files...",SID_LOG_OPEN|SID_LOG_TIMER);
     int        i_g_last;
     halo_properties_SAGE_info *tree;
     tree=(halo_properties_SAGE_info *)SID_calloc(sizeof(halo_properties_SAGE_info)*n_halos_biggest_tree);
     halo=tree;
     i_report   =1;
     next_report=(int)((float)i_report*0.1*(float)n_trees_total);
     for(i_file=0,i_g_last=-1,fp_out=NULL,k_tree=0;i_file<n_in;i_file++){
        int  n_trees_file;
        int  n_halos_file;
        int  i_tree;
        int *n_halos_tree;
        // Open input file
        sprintf(filename_in, "%s/vertical/%sgroup_trees_%03d.dat",filename_root_in,group_text_prefix,i_file);
        fp_in=fopen(filename_in,"r");
        // Read header
        fread(&n_trees_file,sizeof(int),1,fp_in);
        fread(&n_halos_file,sizeof(int),1,fp_in);
        n_halos_tree=(int *)SID_malloc(sizeof(int)*n_trees_file);
        fread(n_halos_tree,sizeof(int),n_trees_file,fp_in);
        for(i_tree=0;i_tree<n_trees_file;i_tree++,k_tree++){
           if(n_halos_tree[i_tree]>0){
              int i_x,i_y,i_z,i_g;
              // Read this tree
              fread(tree,sizeof(halo_properties_SAGE_info),n_halos_tree[i_tree],fp_in);
              // Decide which grid cell it belongs to
              i_x=(int)((float)grid_size*(halo->pos[0]/box_size));i_x=MIN(i_x,grid_size-1);
              i_y=(int)((float)grid_size*(halo->pos[1]/box_size));i_y=MIN(i_y,grid_size-1);
              i_z=(int)((float)grid_size*(halo->pos[2]/box_size));i_z=MIN(i_z,grid_size-1);
              i_g=(i_z*grid_size+i_y)*grid_size+i_x;
              // Perform a sanity check
              if(n_halos_tree[i_tree]!=tree_size[i_g][n_trees_file_out[i_g]])
                 SID_trap_error("Tree size doen't make sense (ie. %d!=%d)",ERROR_LOGIC,n_halos_tree[i_tree],tree_size[i_g][n_trees_file_out[i_g]++]);
              n_trees_file_out[i_g]++;
              // Open the output file
              if(i_g_last!=i_g){
                 sprintf(filename_out,"%s/vertical_grid/%sgroup_trees_%03d.dat",filename_root_in,group_text_prefix,i_g);
                 i_g_last=i_g;
                 if(fp_out!=NULL)
                   fclose(fp_out);
                 fp_out=fopen(filename_out,"a");
              }
              // Write this tree
              fwrite(tree,sizeof(halo_properties_SAGE_info),n_halos_tree[i_tree],fp_out);
              // Status message
              if(k_tree==next_report){
                 SID_log("%3d%% complete.",SID_LOG_COMMENT|SID_LOG_TIMER,i_report*10);
                 i_report++;
                 if(i_report==10)
                    next_report=n_trees_total-1;
                 else
                    next_report=(int)((float)i_report*0.1*(float)n_trees_total);
              }
           }
        }
        fclose(fp_in);
        SID_free(SID_FARG n_halos_tree);
     }
     if(fp_out!=NULL)
        fclose(fp_out);
     SID_free(SID_FARG tree);
     SID_log("Done.",SID_LOG_CLOSE);

     // Clean-up
     for(i_file=0;i_file<n_files_out;i_file++)
        SID_free(SID_FARG tree_size[i_file]);
     SID_free(SID_FARG tree_size);
     SID_free(SID_FARG n_halos_file_out);
     SID_free(SID_FARG n_trees_file_out);

     SID_log("Done.",SID_LOG_CLOSE);
  
  }
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

