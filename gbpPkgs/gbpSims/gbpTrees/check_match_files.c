#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

int main(int argc, char *argv[]){
  int     n_search;
  int     i_halo;
  char    filename_in[MAX_FILENAME_LENGTH];
  char    group_text_prefix[4];
  int     n_files;
  int     k_read;
  int     l_read;
  int    *n_particles_i;
  int    *n_particles_j;
  int     j_read;
  int     mode;
  int     j_halo;
  int     i_read;
  int     i_read_start;
  int     i_read_stop;
  SID_fp  fp_in;

  SID_init(&argc,&argv,NULL,NULL);

  // Fetch user inputs
  char filename_SSimPL_root[MAX_FILENAME_LENGTH];
  strcpy(filename_SSimPL_root,argv[1]);
  SID_log("Checking the integrity of the match files for {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_SSimPL_root);
  int *n_groups   =NULL;
  int *n_subgroups=NULL;
  for(int i_type=0;i_type<2;i_type++){
     // Convert filename_root to filename
     switch(i_type){
        case 0:
        mode=MATCH_SUBGROUPS;
        sprintf(group_text_prefix,"sub");
        break;
        case 1:
        mode=MATCH_GROUPS;
        sprintf(group_text_prefix,"");
        break;
     }
     SID_log("Processing %sgroups...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);

     // Set the standard SSiMPL match file path
     char filename_root_in[MAX_FILENAME_LENGTH];
     sprintf(filename_root_in,"%s/trees/matches/",filename_SSimPL_root);

     // Read halo sizes from header file
     SID_log("Processing header file...",SID_LOG_OPEN|SID_LOG_TIMER);
     sprintf(filename_in,"%s/%sgroup_matches_header.dat",filename_root_in,group_text_prefix);
     SID_fopen(filename_in,"r",&fp_in);
     SID_fread(&i_read_start,sizeof(int),1,&fp_in);
     SID_fread(&i_read_stop, sizeof(int),1,&fp_in);
     SID_fread(&n_search,    sizeof(int),1,&fp_in);
     SID_fread(&n_files,     sizeof(int),1,&fp_in);
     int *n_halos=NULL;
     switch(mode){
        case MATCH_SUBGROUPS:
           n_subgroups=(int *)SID_malloc(sizeof(int)*n_files);
           n_halos    =n_subgroups;
           break;
        case MATCH_GROUPS:
           n_groups=(int *)SID_malloc(sizeof(int)*n_files);
           n_halos =n_groups;
           break;
     }
     if(mode==MATCH_GROUPS) SID_log("Halo counts (snap/No. groups/No. subgroups):",SID_LOG_OPEN);
     for(k_read=0;k_read<n_files;k_read++){
        SID_fread(&l_read,           sizeof(int),1,              &fp_in);
        SID_fread(&(n_halos[k_read]),sizeof(int),1,              &fp_in);
        SID_fskip(                   sizeof(int),n_halos[k_read],&fp_in);
        if(mode==MATCH_GROUPS){
           int *n_subgroups_group=(int *)SID_malloc(sizeof(int)*n_halos[k_read]);
           SID_fread(n_subgroups_group,sizeof(int),n_halos[k_read],&fp_in);
           int n_subgroups_test=0;
           for(int i_test=0;i_test<n_halos[k_read];i_test++)
              n_subgroups_test+=n_subgroups_group[i_test];
           if(n_subgroups[k_read]!=n_subgroups_test)
              SID_log("Error in %s header: l_read=%3d k_read=%3d n_subgroups: %d!=%d\n",SID_LOG_COMMENT,l_read,k_read,n_subgroups[k_read],n_subgroups_test);
           SID_free(SID_FARG n_subgroups_group);
        }
        if(mode==MATCH_GROUPS) SID_log("%03d %d %d",SID_LOG_COMMENT,k_read,n_groups[k_read],n_subgroups[k_read]);
     }
     if(mode==MATCH_GROUPS) SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);
     SID_fclose(&fp_in);
     SID_log("Done.",SID_LOG_CLOSE);

     SID_log("Processing match files...",SID_LOG_OPEN|SID_LOG_TIMER);
     for(int i_read=i_read_start;i_read<i_read_stop;i_read++){
        for(int j_read=MAX(0,i_read-n_search);j_read<MIN(i_read_stop,i_read+n_search);j_read++){
           if(i_read!=j_read){
              sprintf(filename_in,"%s/%03d/%sgroup_matches_%03d_%03d.dat",filename_root_in,i_read,group_text_prefix,i_read,j_read);
              SID_log("Processing {%s}...",SID_LOG_OPEN,filename_in);

              // Read header information
              int i_read_in;
              int j_read_in;
              int n_groups_i;
              int n_groups_j;
              SID_fopen(filename_in,"r",&fp_in);
              SID_fread(&i_read_in, sizeof(int),1,&fp_in);
              SID_fread(&j_read_in, sizeof(int),1,&fp_in);
              SID_fread(&n_groups_i,sizeof(int),1,&fp_in);
              SID_fread(&n_groups_j,sizeof(int),1,&fp_in);

              if(i_read_in!=i_read || j_read_in!=j_read || n_groups_i!=n_halos[n_files-i_read_in-1] || n_groups_j!=n_halos[n_files-j_read_in-1])
                 SID_log("Error in matching file: i_read=%3d j_read=%3d n_i_in=%d n_i=%d n_j_in=%d n_j=%d\n",SID_LOG_COMMENT,i_read,j_read,n_groups_i,n_halos[n_files-i_read_in-1],n_groups_j,n_halos[n_files-j_read_in-1]);
         
              // Read matches
              int match;
              for(k_read=0;k_read<n_groups_i;k_read++)
                 SID_fread(&match,sizeof(int),1,&fp_in);
         
              // Read indices
              size_t indices;
              for(k_read=0;k_read<n_groups_i;k_read++)
                 SID_fread(&indices,sizeof(size_t),1,&fp_in);
         
              // Read scores
              float score;
              for(k_read=0;k_read<n_groups_i;k_read++)
                 SID_fread(&score,sizeof(float),1,&fp_in);
         
              // Close file
              SID_fclose(&fp_in);

              SID_log("Done.",SID_LOG_CLOSE);
           }
        }
     }
     SID_log("Done.",SID_LOG_CLOSE);

     SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_free(SID_FARG n_groups);
  SID_free(SID_FARG n_subgroups);

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

