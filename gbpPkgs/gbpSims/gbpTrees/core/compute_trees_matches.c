#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>

void check_for_tree_matches_local(char *filename_root_out,int i_read_start,int i_read_stop,int i_read,int n_search_total,int *flag_go,int **flag_go_array);
void check_for_tree_matches_local(char *filename_root_out,int i_read_start,int i_read_stop,int i_read,int n_search_total,int *flag_go,int **flag_go_array){
  int   j_read;
  int   k_read;
  int   k_match;
  int   k_order;
  int   i_read_order;
  int   j_read_order;
  char  filename_cat1_order[256];
  char  filename_cat2_order[256];
  char  group_text_prefix[5];
  char  filename_out[256];
  char  filename_out_dir[256];
  char  filename_out_name[256];
  char  filename_out_dir_snap[256];
  FILE *fp_test;

  if((*flag_go_array)==NULL)
     (*flag_go_array)=(int *)SID_malloc(sizeof(int)*n_search_total);
  for(k_read=0;k_read<n_search_total;k_read++)
     (*flag_go_array)[k_read]=FALSE;

  if(SID.I_am_Master){
     strcpy(filename_out_dir, filename_root_out);
     strcpy(filename_out_name,filename_root_out);
     strip_path(filename_out_name);
     for(j_read=i_read-1,k_read=0;j_read>=i_read-n_search_total;j_read--,k_read++){
        if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
           for(k_order=0;k_order<2;k_order++){
              if(k_order==0){
                 i_read_order=i_read;
                 j_read_order=j_read;
              }
              else{
                 i_read_order=j_read;
                 j_read_order=i_read;
              }
              sprintf(filename_cat1_order,"%03d",i_read_order);
              sprintf(filename_cat2_order,"%03d",j_read_order);
              sprintf(filename_out_dir_snap,"%s/%s",filename_out_dir,filename_cat1_order);
              for(k_match=0;k_match<2;k_match++){
                 switch(k_match){
                    case 0:
                    sprintf(group_text_prefix,"sub");
                    break;
                    case 1:
                    sprintf(group_text_prefix,"");
                    break;
                 }
                 if(filename_out_dir!=NULL)
                    sprintf(filename_out,"%s/%sgroup_matches_%s_%s.dat",filename_out_dir_snap,group_text_prefix,filename_cat1_order,filename_cat2_order);
                 else
                    sprintf(filename_out,"%s_%sgroup_matches_%s_%s.dat",filename_out_name,    group_text_prefix,filename_cat1_order,filename_cat2_order);
                 if((fp_test=fopen(filename_out,"r"))!=NULL)
                    fclose(fp_test);
                 else{
                    (*flag_go)              =TRUE;
                    (*flag_go_array)[k_read]=TRUE;
                 }
              }
           }
        }
     }
  }
  SID_Bcast(flag_go,         sizeof(int),               MASTER_RANK,SID.COMM_WORLD);  
  SID_Bcast((*flag_go_array),sizeof(int)*n_search_total,MASTER_RANK,SID.COMM_WORLD);  
}

int compute_trees_matches(char   *filename_root_in,
                          char   *filename_root_out,
                          int     i_read_start,
                          int     i_read_stop,
                          int     i_read_step,
                          int     n_search,
                          int     mode){
  char        filename_out[256];
  char        filename_out_dir[256];
  char        filename_out_name[256];
  FILE       *fp_test;
  FILE       *fp_out;
  SID_fp      fp_in;
  int         i_read,k_read,l_read;
  int         flag_go;
  int         i_read_start_file;
  int         i_read_stop_file;
  int         i_read_step_file;
  int         n_search_file;
  int         k_match;
  int         n_search_total;
  int         n_k_match;
  int         flag_match_subgroups;
  char        group_text_prefix[5];
  int         n_matches;
  int         n_files;
  int         j_read;
  char        filename_cat1[256];
  char        filename_cat2[256];
  char        filename_cat1_order[256];
  char        filename_cat2_order[256];
  char        filename_out_dir_snap[256];
  int         n_groups_1;
  int         n_groups_1_local;
  int         n_groups_2;
  int         n_groups_2_local;
  int         i_group;
  int         buffered_count;
  int         buffered_count_local;
  int         j_group;
  int         index_test;
  int         i_rank;
  int        *n_particles;
  int        *n_sub_group;
  int        *match_id   =NULL;
  float      *match_score=NULL;
  char        cat_name_1[20];
  char        cat_name_2[20];
  size_t     *match_rank =NULL;
  size_t     *match_index=NULL;
  size_t      offset;
  plist_info  plist1;
  plist_info  plist2;
  void       *buffer;
  int        *buffer_int;
  size_t     *buffer_size_t;
  float      *buffer_float;
  int         n_buffer_max=1000;
  int         n_buffer;
  int         i_buffer;
  int         j_buffer;
  int         flag_compute_header_subgroups;
  int         flag_compute_header_groups;
  int         flag_compute_header;
  int         flag_sucessful_completion=TRUE;

  if(check_mode_for_flag(mode,WRITE_MATCHES_PERFORM_CHECK))
     SID_log("Validating merger tree matches...",SID_LOG_OPEN|SID_LOG_TIMER);
  else
     SID_log("Constructing merger tree matches...",SID_LOG_OPEN|SID_LOG_TIMER);

  n_search_total=n_search*i_read_step;

  // Set results directory and root
  strcpy(filename_out_dir, filename_root_out);
  strcpy(filename_out_name,filename_root_out);
  strip_path(filename_out_name);

  // Check if all the needed files are present
  int flag_all_inputs_present=TRUE;
  for(i_read=i_read_stop;i_read>=i_read_start && flag_all_inputs_present;i_read--)
     flag_all_inputs_present&=check_for_matching_input_files(filename_root_in,i_read);

  // Generate headers if needed and possible
  int flag_create_headers;
  if(flag_all_inputs_present)
     flag_create_headers=check_mode_for_flag(mode,WRITE_MATCHES_CHECK_HEADER) ||
                         check_mode_for_flag(mode,WRITE_MATCHES_PERFORM_CHECK);
  else{
     flag_create_headers      =FALSE;
     flag_sucessful_completion=FALSE;
  }
  if(flag_create_headers){
     SID_log("Checking if header files exist and are adequate...",SID_LOG_OPEN|SID_LOG_TIMER);
     for(k_match=0;k_match<2;k_match++){

        // Set working array pointers to point to group or subgroup information (alternately, depending on k_match)
        switch(k_match){
           case 0:
           flag_match_subgroups=MATCH_SUBGROUPS;
           sprintf(group_text_prefix,"sub");
           break;
           case 1:
           flag_match_subgroups=MATCH_GROUPS;
           sprintf(group_text_prefix,"");
           break;
        }

        // Set filenames
        sprintf(filename_out,"%s/%sgroup_matches_header.dat",filename_out_dir,group_text_prefix);
        flag_go=TRUE;
        if((fp_test=fopen(filename_out,"r"))!=NULL){
           fclose(fp_test);
           SID_fopen(filename_out,"r",&fp_in);
           SID_fread_all(&i_read_start_file,sizeof(int),1,&fp_in);
           SID_fread_all(&i_read_stop_file, sizeof(int),1,&fp_in);
           SID_fread_all(&n_search_file,    sizeof(int),1,&fp_in);
           SID_fread_all(&n_files,          sizeof(int),1,&fp_in);
           SID_fclose(&fp_in);
           if(i_read_stop_file<i_read_stop || i_read_start_file>i_read_start || n_search_file<n_search_total){
              if(k_match==0)
                 flag_compute_header_subgroups=TRUE;
              else
                 flag_compute_header_groups=TRUE;
              SID_log("Header file for %sgroups needs to be computed.",SID_LOG_COMMENT,group_text_prefix);
           }
           else{
              if(k_match==0)
                 flag_compute_header_subgroups=FALSE;
              else
                 flag_compute_header_groups=FALSE;
              SID_log("Header file for %sgroups is fine.",SID_LOG_COMMENT,group_text_prefix);
           }
        }
        else{
           if(k_match==0)
              flag_compute_header_subgroups=TRUE;
           else
              flag_compute_header_groups   =TRUE;
        }
     }
     SID_log("Done.",SID_LOG_CLOSE);

     // Write the header ...
     if(flag_compute_header_subgroups || flag_compute_header_groups){
        SID_log("Constructing header files...",SID_LOG_OPEN|SID_LOG_TIMER);
        SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

        // ... count the number of matches and files that we will store ...
        for(i_read=i_read_stop,n_matches=0,n_files=0;i_read>=i_read_start;i_read--){
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 n_matches++;
              }
           }
           n_files++;
        }

        // ... loop over each snapshot ...
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           sprintf(filename_cat1,"%03d",i_read);
           init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
           read_groups(filename_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_NOIDS,&plist1,filename_cat1);
           for(k_match=0;k_match<2;k_match++){
              switch(k_match){
                 case 0:
                 flag_match_subgroups=MATCH_SUBGROUPS;
                 sprintf(group_text_prefix,"sub");
                 flag_compute_header=flag_compute_header_subgroups;
                 break;
                 case 1:
                 flag_match_subgroups=MATCH_GROUPS;
                 sprintf(group_text_prefix,"");
                 flag_compute_header=flag_compute_header_groups;
                 break;
              }
              if(flag_compute_header){
                 sprintf(filename_out,"%s/%sgroup_matches_header.dat",filename_out_dir,group_text_prefix);

                 // This gets written just once 
                 if(i_read==i_read_stop && SID.I_am_Master){
                    mkdir(filename_out_dir,02755);
                    fp_out=fopen(filename_out,"w");
                    fwrite(&i_read_start,  sizeof(int),1,fp_out);
                    fwrite(&i_read_stop,   sizeof(int),1,fp_out);
                    fwrite(&n_search_total,sizeof(int),1,fp_out);
                    fwrite(&n_files,       sizeof(int),1,fp_out);
                    fclose(fp_out);
                 }

                 // This gets written for every snapshot
                 n_groups_1      =((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_all_%s",group_text_prefix,filename_cat1))[0];
                 n_groups_1_local=((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_%s",    group_text_prefix,filename_cat1))[0];
                 if(SID.I_am_Master){
                    fp_out=fopen(filename_out,"a");
                    fwrite(&i_read,    sizeof(int),1,fp_out);
                    fwrite(&n_groups_1,sizeof(int),1,fp_out);
                    fclose(fp_out);
                 }
                 SID_Barrier(SID.COMM_WORLD);

                 // Write the arrays for each snapshot
                 fp_out=fopen(filename_out,"a");
                 if(n_groups_1>0){

                    // Initialize a buffer for writing
                    int   n_groups_1_local_max;
                    int   n_buffer;
                    char *buffer;
                    calc_max_global(&n_groups_1_local,&n_groups_1_local_max,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
                    buffer=(char *)SID_malloc(n_groups_1_local_max*sizeof(int));

                    // Write the group/subgroup sizes
                    int i_rank;
                    n_particles=(int *)ADaPS_fetch(plist1.data,"n_particles_%sgroup_%s",group_text_prefix,filename_cat1);
                    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
                       if(i_rank==SID.My_rank){
                          n_buffer=n_groups_1_local;
                          memcpy(buffer,n_particles,n_buffer*sizeof(int));
                       }
                       SID_Bcast(&n_buffer,         sizeof(int),i_rank,SID.COMM_WORLD);
                       SID_Bcast(buffer,   n_buffer*sizeof(int),i_rank,SID.COMM_WORLD);
                       if(SID.I_am_Master)
                          fwrite(buffer,n_buffer,sizeof(int),fp_out);
                    }

                    // Write the number of substructures per group
                    if(flag_match_subgroups==MATCH_GROUPS){
                       n_sub_group=(int *)ADaPS_fetch(plist1.data,"n_subgroups_group_%s",filename_cat1);
                       for(i_rank=0;i_rank<SID.n_proc;i_rank++){
                          if(i_rank==SID.My_rank){
                             n_buffer=n_groups_1_local;
                             memcpy(buffer,n_sub_group,n_buffer*sizeof(int));
                          }
                          SID_Bcast(&n_buffer,         sizeof(int),i_rank,SID.COMM_WORLD);
                          SID_Bcast(buffer,   n_buffer*sizeof(int),i_rank,SID.COMM_WORLD);
                          if(SID.I_am_Master)
                             fwrite(buffer,n_buffer,sizeof(int),fp_out);
                       }
                    }

                    // Free the buffer
                    SID_free(SID_FARG buffer);
                 }
                 fclose(fp_out);
              } // If flag_compute_header
           } // Loop over k_match
           free_plist(&plist1);
        }
        SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
        SID_log("Done.",SID_LOG_CLOSE);
     }
  }

  // Check to see if there are any matches needing to be completed
  i_read =i_read_stop;
  //flag_go=TRUE; // Uncomment this to force recalulation of header files
  flag_go=FALSE;
  int *flag_go_array=NULL;
  SID_log("Checking for matching files...",SID_LOG_OPEN);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  for(;i_read>i_read_start && !flag_go;i_read--)
     check_for_tree_matches_local(filename_out_dir,i_read_start,i_read_stop,i_read,n_search_total,&flag_go,&flag_go_array);
  if(flag_go)
     i_read++;
  if(!flag_go)
     SID_log("All matches present and complete.",SID_LOG_COMMENT);
  else if(i_read!=i_read_stop)
     SID_log("Matching only needs to start at snapshot #%03d.",SID_LOG_COMMENT,i_read);
  else
     SID_log("Matching starting with first snapshot.",SID_LOG_COMMENT);
  SID_log("Done.",SID_LOG_CLOSE);

  // Generate matches (if needed).  Loop over base groups first...
  if(flag_go){
     if(check_mode_for_flag(mode,WRITE_MATCHES_PERFORM_CHECK))
        SID_log("Checking that all needed matches are present...",SID_LOG_OPEN|SID_LOG_TIMER);
     else
        SID_log("Generating matches...",SID_LOG_OPEN|SID_LOG_TIMER);
     SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,1);
     // i_read gets initialized by the previous loop
     for(;i_read>i_read_start;i_read--){
        flag_go=FALSE;
        check_for_tree_matches_local(filename_out_dir,i_read_start,i_read_stop,i_read,n_search_total,&flag_go,&flag_go_array);
        if(flag_go){
           if(check_for_matching_input_files(filename_root_in,i_read)){
              SID_log("Processing matches for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);                 

              // Read base group
              sprintf(filename_cat1,"%03d",i_read);
              init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
              SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
              read_groups(filename_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,&plist1,filename_cat1,-1,-1);
              SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

              // Compute each matching combaination and write the results to the file
              int         k_read;
              int         i_read_order;
              int         j_read_order;
              plist_info *plist1_order;
              plist_info *plist2_order;
              for(j_read=i_read-1,k_read=0;j_read>=i_read-n_search_total;j_read--,k_read++){
                 if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                    if(check_for_matching_input_files(filename_root_in,j_read)){
                       int k_order;
                       int flag_read;
                       // Loop over forward/back matching
                       SID_log("Processing matches between snaps %03d and %03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read,j_read);
                       sprintf(filename_cat2,"%03d",j_read);

                       // If this snapshot combination is missing at least one file, proceed.
                       flag_read=TRUE;
                       if(flag_go_array[k_read]){
                          for(k_order=0;k_order<2;k_order++){
                             if(k_order==0){
                                i_read_order=i_read;
                                j_read_order=j_read;
                                plist1_order=&plist1;
                                plist2_order=&plist2;
                             }
                             else{
                                i_read_order=j_read;
                                j_read_order=i_read;
                                plist1_order=&plist2;
                                plist2_order=&plist1;
                             }
                             sprintf(filename_cat1_order,"%03d",i_read_order);
                             sprintf(filename_cat2_order,"%03d",j_read_order);

                             // Read catalog to match to
                             if(flag_read){
                                int PHK_min_local;
                                int PHK_max_local;
                                PHK_min_local=((int *)ADaPS_fetch(plist1.data,"PHK_min_local_%s",filename_cat1))[0];
                                PHK_max_local=((int *)ADaPS_fetch(plist1.data,"PHK_max_local_%s",filename_cat1))[0];
                                init_plist(&plist2,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
                                SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
                                read_groups(filename_root_in,j_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,
                                            &plist2,filename_cat2,PHK_min_local,PHK_max_local);
                                SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
                                flag_read=FALSE;
                             }

                             // Perform matching
                             for(k_match=0;k_match<2;k_match++){
                                switch(k_match){
                                   case 0:
                                   flag_match_subgroups=MATCH_SUBGROUPS;
                                   sprintf(group_text_prefix,"sub");
                                   break;
                                   case 1:
                                   flag_match_subgroups=MATCH_GROUPS;
                                   sprintf(group_text_prefix,"");
                                   break;
                                }
                                SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
                                match_halos(plist1_order,
                                            NULL,
                                            i_read_order,
                                            NULL,
                                            0,
                                            plist2_order,
                                            NULL,
                                            j_read_order,
                                            NULL,
                                            0,
                                            "match",
                                            flag_match_subgroups|MATCH_STORE_SCORE,
                                            MATCH_SCORE_RANK_INDEX);
                                SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

                                // Writing results
                                SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
                                write_match_results(filename_out_dir,
                                                    filename_out_name,
                                                    i_read_order,
                                                    j_read_order,
                                                    filename_cat1_order,
                                                    filename_cat2_order,
                                                    plist1_order,
                                                    plist2_order,
                                                    k_match,
                                                    WRITE_MATCHES_MODE_TREES);
                                SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
                             }
                          } // Loop over matching order
                       } // If this match result didn't already exist
                       else if(!check_mode_for_flag(mode,WRITE_MATCHES_PERFORM_CHECK))
                          SID_log("Matching for %03d->%03d present...Skipping.",SID_LOG_COMMENT,i_read,j_read);
                       if(!flag_read)
                          free_plist(&plist2);
                       SID_log("Done.",SID_LOG_CLOSE);
                    }
                    else{
                       SID_log("Input for snapshot #%d {%s} absent...Skipping.",SID_LOG_COMMENT,j_read,filename_root_in);
                       flag_sucessful_completion=FALSE;
                    }
                 } // If this is a valid pair
              } // Second loop over snapshots
              free_plist(&plist1);
              SID_log("Done.",SID_LOG_CLOSE);
           }
           else{
              SID_log("Input for snapshot #%d {%s} absent...Skipping.",SID_LOG_COMMENT,i_read,filename_root_in);                 
              flag_sucessful_completion=FALSE;
           } 
        } // If this snapshot needs to be processed
        else if(!check_mode_for_flag(mode,WRITE_MATCHES_PERFORM_CHECK))
           SID_log("Matching for snapshot #%d present...Skipping.",SID_LOG_COMMENT,i_read);
     } // Loop over base snapshots
     SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Clean-up
  SID_free(SID_FARG flag_go_array);

  SID_log("Done.",SID_LOG_CLOSE);
  return(flag_sucessful_completion);
}

