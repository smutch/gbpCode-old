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

void compute_trees_matches(char   *filename_root_in,
                           char   *filename_root_out,
                           int     i_read_start,
                           int     i_read_stop,
                           int     i_read_step,
                           int    *n_files_return,
                           int   **n_subgroups_return,
                           int   **n_groups_return,
                           int     n_search){
  char        filename_out[256];
  char        filename_out_dir[256];
  char        filename_out_name[256];
  FILE       *fp_test;
  SID_fp      fp_out;
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
  int        *n_return;
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

  SID_log("Constructing merger tree matches...",SID_LOG_OPEN|SID_LOG_TIMER);

  n_search_total=n_search*i_read_step;

  // Initialize the n_halos arrays that are to be returned
  for(i_read=i_read_stop,(*n_files_return)=0;      
      i_read>=i_read_start;
      i_read-=i_read_step) (*n_files_return)++;  
  (*n_subgroups_return)=(int *)SID_malloc(sizeof(int)*(*n_files_return));
  (*n_groups_return)   =(int *)SID_malloc(sizeof(int)*(*n_files_return));

  // Set results directory and root
  strcpy(filename_out_dir, filename_root_out);
  strcpy(filename_out_name,filename_root_out);
  strip_path(filename_out_name);

  // Generate headers if needed
  int flag_create_headers=FALSE;
  if(flag_create_headers){
     SID_log("Checking if header files exist and are adequate...",SID_LOG_OPEN);
     for(k_match=0;k_match<2;k_match++){

        // Set working array pointers to point to group or subgroup information (alternately, depending on k_match)
        switch(k_match){
           case 0:
           flag_match_subgroups=MATCH_SUBGROUPS;
           sprintf(group_text_prefix,"sub");
           n_return=(*n_subgroups_return);
           break;
           case 1:
           flag_match_subgroups=MATCH_GROUPS;
           sprintf(group_text_prefix,"");
           n_return=(*n_groups_return);
           break;
        }

        // Set filenames
        sprintf(filename_out,"%s/%s.%sgroup_matches_header",filename_out_dir,filename_out_name,group_text_prefix);

        flag_go=TRUE;
        if((fp_test=fopen(filename_out,"r"))!=NULL){
           fclose(fp_test);
           SID_fopen(filename_out,"r",&fp_in);
           SID_fread_all(&i_read_start_file,sizeof(int),1,&fp_in);
           SID_fread_all(&i_read_stop_file, sizeof(int),1,&fp_in);
           SID_fread_all(&n_search_file,    sizeof(int),1,&fp_in);
           SID_fread_all(&n_files,          sizeof(int),1,&fp_in);
           SID_fclose(&fp_in);
           if(i_read_stop_file>=i_read_stop && i_read_start_file<=i_read_start && n_search_file>=n_search_total){
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
     if(flag_compute_header_subgroups || flag_compute_header_groups)
        SID_log("Constructing header files...",SID_LOG_OPEN|SID_LOG_TIMER);
     //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
     SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

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
              n_return=(*n_subgroups_return);
              flag_compute_header=flag_compute_header_subgroups;
              break;
              case 1:
              flag_match_subgroups=MATCH_GROUPS;
              sprintf(group_text_prefix,"");
              n_return=(*n_groups_return);
              flag_compute_header=flag_compute_header_groups;
              break;
           }
           if(flag_compute_header){
              sprintf(filename_out,"%s/%s.%sgroup_matches_header",filename_out_dir,filename_out_name,group_text_prefix);
      
              if(i_read==i_read_stop && SID.I_am_Master){
                 mkdir(filename_out_dir,02755);
                 SID_fopen(filename_out,"w",&fp_out);
                 SID_fwrite(&i_read_start,  sizeof(int),1,&fp_out);
                 SID_fwrite(&i_read_stop,   sizeof(int),1,&fp_out);
                 SID_fwrite(&n_search_total,sizeof(int),1,&fp_out);
                 SID_fwrite(&n_files,       sizeof(int),1,&fp_out);
                 SID_fclose(&fp_out);
              }

              n_groups_1      =((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_all_%s",group_text_prefix,filename_cat1))[0];
              n_groups_1_local=((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_%s",    group_text_prefix,filename_cat1))[0];
              if(SID.I_am_Master){
                 SID_fopen(filename_out,"a",&fp_out);
                 SID_fwrite(&i_read,    sizeof(int),1,&fp_out);
                 SID_fwrite(&n_groups_1,sizeof(int),1,&fp_out);
                 SID_fclose(&fp_out);
              }
              SID_Barrier(SID.COMM_WORLD);
              SID_fopen(filename_out,"a",&fp_out);
              if(n_groups_1>0){
                 n_particles=(int *)ADaPS_fetch(plist1.data,"n_particles_%sgroup_%s",group_text_prefix,filename_cat1);
                 SID_fwrite_ordered(n_particles,sizeof(int),(size_t)n_groups_1_local,&fp_out);
                 if(flag_match_subgroups==MATCH_GROUPS){
                    n_sub_group=(int *)ADaPS_fetch(plist1.data,"n_subgroups_group_%s",filename_cat1);
                    SID_fwrite_ordered(n_sub_group,sizeof(int),(size_t)n_groups_1_local,&fp_out);
                 }
              }
              SID_fclose(&fp_out);
           } // If flag_compute_header
        } // Loop over k_match
        free_plist(&plist1);
     }
     SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
     if(flag_compute_header_subgroups || flag_compute_header_groups)
        SID_log("Done.",SID_LOG_CLOSE);
  }

  // Generate matches.  Loop over base groups first...
  SID_log("Generating matches...",SID_LOG_OPEN);
  SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  for(i_read=i_read_stop;i_read>i_read_start;i_read--){
     SID_log("Processing matches for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);                 

     // Read base group
     sprintf(filename_cat1,"%03d",i_read);
     init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
     SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
     read_groups(filename_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,&plist1,filename_cat1,-1,-1);
     SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

     // Compute each matching combaination and write the results to the file
     int         i_read_order;
     int         j_read_order;
     plist_info *plist1_order;
     plist_info *plist2_order;
     for(j_read=i_read-1;j_read>=i_read-n_search_total;j_read--){
        if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
           int k_order;
           int flag_read;
           // Loop over forward/back matching
           SID_log("Processing matches between snaps %03d and %03d...",SID_LOG_OPEN,i_read,j_read);                 
           sprintf(filename_cat2,"%03d",j_read);

           // Check for valid files.  If they exist and are ok, then
           //   we don't need to read or process this snapshot
           if(SID.I_am_Master){
              flag_go=FALSE;
              for(k_order=0,flag_read=TRUE;k_order<2;k_order++){
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
                 for(k_match=0;k_match<2;k_match++){
                    switch(k_match){
                       case 0:
                       sprintf(group_text_prefix,"sub");
                       break;
                       case 1:
                       sprintf(group_text_prefix,"");
                       break;
                    }
                    sprintf(filename_out,"%s/%s_%s_%s.%sgroup_matches",
                                         filename_out_dir,
                                         filename_out_name,
                                         filename_cat1_order,
                                         filename_cat2_order,
                                         group_text_prefix);
                    if((fp_test=fopen(filename_out,"r"))!=NULL)
                       fclose(fp_test);
                    else
                       flag_go=TRUE;
                 }
              }
           }
           SID_Bcast(&flag_go,sizeof(int),MASTER_RANK,SID.COMM_WORLD);

           // If this snapshot combination is missing at least one file, proceed.
           if(flag_go){
              for(k_order=0,flag_read=TRUE;k_order<2;k_order++){
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
                    match_halos(plist1_order,i_read_order,NULL,0,plist2_order,j_read_order,NULL,0,"match",flag_match_subgroups|MATCH_STORE_SCORE);

                    // Writing results
                    write_match_results(filename_out_dir,filename_out_name,i_read_order,j_read_order,plist1_order,plist2_order,k_match);
                 }
              } // Loop over matching order
           } // If this match result didn't already exist
           else
              SID_log("Skipping %03d->%03d matching.",SID_LOG_COMMENT,i_read_order,j_read_order);                 
           if(!flag_read)
              free_plist(&plist2);
           SID_log("Done.",SID_LOG_CLOSE);
        } // If this is a valid pair
     } // Second loop over snapshots
     free_plist(&plist1);
     SID_log("Done.",SID_LOG_CLOSE);
  } // Loop over base snapshots
  SID_log("Done.",SID_LOG_CLOSE);

     // Construct the n_halos arrays that are to be returned
     /*
     SID_fopen(filename_out,"r",&fp_in);
     SID_fread_all(&i_read_start_file,sizeof(int),1,&fp_in);
     SID_fread_all(&i_read_stop_file, sizeof(int),1,&fp_in);
     SID_fread_all(&n_search,         sizeof(int),1,&fp_in);
     SID_fread_all(&n_files,          sizeof(int),1,&fp_in);
     for(i_read=i_read_stop_file,j_read=i_read,k_read=0;i_read>=i_read_start_file;i_read--){
        SID_fread_all(&l_read,    sizeof(int),1,&fp_in);
        SID_fread_all(&n_groups_1,sizeof(int),1,&fp_in);
        if(flag_match_subgroups==MATCH_GROUPS)
           SID_fseek(&fp_in,2*sizeof(int),n_groups_1,SID_SEEK_CUR);
        else
           SID_fseek(&fp_in,sizeof(int),n_groups_1,SID_SEEK_CUR);
        if(j_read==l_read){
           n_return[k_read++]=n_groups_1;
           j_read-=i_read_step;
        }
     }
     SID_fclose(&fp_in);
    */
  SID_log("Done.",SID_LOG_CLOSE);
}

