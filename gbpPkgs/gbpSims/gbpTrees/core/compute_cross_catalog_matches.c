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

int compute_cross_catalog_matches(char   *filename_root_in_1,
                                  char   *filename_root_in_2,
                                  char   *filename_root_out,
                                  int     i_read_1,
                                  int     i_read_2,
                                  float   match_weight_rank_index){
  char        filename_out[256];
  char       *filename_out_dir;
  char        filename_out_name[256];
  FILE       *fp_test;
  FILE       *fp_out;
  SID_fp      fp_in;
  int         k_read,l_read;
  int         flag_go;
  int        *flag_go_array=NULL;
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
  char        filename_cat1[256];
  char        filename_cat2[256];
  char        filename_cat1_order[256];
  char        filename_cat2_order[256];
  char       *filename_out_dir_snap;
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

  if(!strcmp(filename_root_in_1,filename_root_in_2))
     SID_log("Processing matches for catalog {%s} bewtween snapshots #%03d and #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,
             filename_root_in_1,
             i_read_1,
             i_read_2);
  else
     SID_log("Processing matches between catalog;snap's {%s;%03d} and {%s;%03d}...",SID_LOG_OPEN|SID_LOG_TIMER,
             filename_root_in_1,
             i_read_1,
             filename_root_in_2,
             i_read_2);

  // Set results directory and root
  filename_out_dir=NULL;
  strcpy(filename_out_name,filename_root_out);

  // Check if all the needed files are present
  int flag_all_inputs_present=TRUE;
  flag_all_inputs_present&=check_for_matching_input_files(filename_root_in_1,i_read_1);
  flag_all_inputs_present&=check_for_matching_input_files(filename_root_in_2,i_read_2);

  // Read groups
  int PHK_min_local;
  int PHK_max_local;
  int n_bits_PHK_1;
  int n_bits_PHK_2;
  if(i_read_1==i_read_2){
     sprintf(filename_cat1,"A");
     sprintf(filename_cat2,"B");
  }
  else{
     sprintf(filename_cat1,"%03d",i_read_1);
     sprintf(filename_cat2,"%03d",i_read_2);
  }
  //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_groups(filename_root_in_1,i_read_1,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,&plist1,filename_cat1,-1,-1);
  PHK_min_local=((int *)ADaPS_fetch(plist1.data,"PHK_min_local_%s",filename_cat1))[0];
  PHK_max_local=((int *)ADaPS_fetch(plist1.data,"PHK_max_local_%s",filename_cat1))[0];
  n_bits_PHK_1 =((int *)ADaPS_fetch(plist1.data,"n_bits_PHK_%s",   filename_cat1))[0];
  init_plist(&plist2,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  read_groups(filename_root_in_2,i_read_2,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,
              &plist2,filename_cat2,PHK_min_local,PHK_max_local);
  n_bits_PHK_2 =((int *)ADaPS_fetch(plist2.data,"n_bits_PHK_%s",   filename_cat2))[0];
  //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  if(n_bits_PHK_1!=n_bits_PHK_2)
     SID_trap_error("The PHK bit-sizeo of the two catalog's precomputed PHKs don't match (ie %d!=%d).",ERROR_LOGIC,n_bits_PHK_1,n_bits_PHK_2);

  // First match one way, then the other
  int         k_order;
  int         i_read_order;
  int         j_read_order;
  plist_info *plist1_order;
  plist_info *plist2_order;
  for(k_order=0;k_order<2;k_order++){
     if(k_order==0){
        i_read_order=i_read_1;
        j_read_order=i_read_2;
        plist1_order=&plist1;
        plist2_order=&plist2;
        if(i_read_order==j_read_order)
           sprintf(filename_out_name,"%s",filename_root_out);
        if(i_read_1==i_read_2){
           sprintf(filename_cat1_order,"A");
           sprintf(filename_cat2_order,"B");
        }
        else{
           sprintf(filename_cat1_order,"%03d",i_read_1);
           sprintf(filename_cat2_order,"%03d",i_read_2);
        }
     }
     else{
        i_read_order=i_read_2;
        j_read_order=i_read_1;
        plist1_order=&plist2;
        plist2_order=&plist1;
        if(i_read_order==j_read_order)
           sprintf(filename_out_name,"%s",filename_root_out);
        if(i_read_1==i_read_2){
           sprintf(filename_cat1_order,"B");
           sprintf(filename_cat2_order,"A");
        }
        else{
           sprintf(filename_cat1_order,"%03d",i_read_2);
           sprintf(filename_cat2_order,"%03d",i_read_1);
        }
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
                    filename_cat1_order,
                    i_read_order,
                    NULL,
                    0,
                    plist2_order,
                    filename_cat2_order,
                    j_read_order,
                    NULL,
                    0,
                    "match",
                    flag_match_subgroups|MATCH_STORE_SCORE,
                    match_weight_rank_index);
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
                            WRITE_MATCHES_MODE_SINGLE);
        SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
     }
  } // Loop over matching order

  // Clean-up
  free_plist(&plist1);
  free_plist(&plist2);

  //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_log("Done.",SID_LOG_CLOSE);
  return(flag_sucessful_completion);
}

