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

void check_for_tree_matches_local(char *filename_root_out,int i_read_start,int i_read_stop,int i_read,int n_search_total,int *flag_go);
void check_for_tree_matches_local(char *filename_root_out,int i_read_start,int i_read_stop,int i_read,int n_search_total,int *flag_go){
  int   j_read;
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

  if(SID.I_am_Master){
     strcpy(filename_out_dir, filename_root_out);
     strcpy(filename_out_name,filename_root_out);
     strip_path(filename_out_name);
     for(j_read=i_read-1;j_read>=i_read-n_search_total && !(*flag_go);j_read--){
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
                    sprintf(filename_out,"%s/%s_%s_%s.%sgroup_matches",filename_out_dir_snap,filename_out_name,filename_cat1_order,filename_cat2_order,group_text_prefix);
                 else
                    sprintf(filename_out,"%s_%s_%s.%sgroup_matches",filename_out_name,filename_cat1_order,filename_cat2_order,group_text_prefix);
                 if((fp_test=fopen(filename_out,"r"))!=NULL)
                    fclose(fp_test);
                 else
                    (*flag_go)=TRUE;
              }
           }
        }
     }
  }
  SID_Bcast(flag_go,sizeof(int),MASTER_RANK,SID.COMM_WORLD);  
}

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

  SID_log("Constructing merger tree matches...",SID_LOG_OPEN|SID_LOG_TIMER);

  n_search_total=n_search*i_read_step;

  // Set results directory and root
  strcpy(filename_out_dir, filename_root_out);
  strcpy(filename_out_name,filename_root_out);
  strip_path(filename_out_name);

  // Generate headers if needed
  int flag_create_headers=TRUE;
  if(flag_create_headers){
     SID_log("Checking if header files exist and are adequate...",SID_LOG_OPEN);
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
        //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

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
                 sprintf(filename_out,"%s/%s.%sgroup_matches_header",filename_out_dir,filename_out_name,group_text_prefix);

                 // Write the header's header ;)         
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
        SID_log("Done.",SID_LOG_CLOSE);
     }
  }

  // Check to see if there are any matches needing to be completed
  SID_log("Checking for matching files...",SID_LOG_OPEN);
  //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
  SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
  flag_go=FALSE;
  for(i_read=i_read_stop;i_read>i_read_start && !flag_go;i_read--)
     check_for_tree_matches_local(filename_out_dir,i_read_start,i_read_stop,i_read,n_search_total,&flag_go);
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
     SID_log("Generating matches...",SID_LOG_OPEN);
     //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
     SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
     // i_read gets initialized by the previous loop
     for(;i_read>i_read_start;i_read--){
        flag_go=FALSE;
        check_for_tree_matches_local(filename_out_dir,i_read_start,i_read_stop,i_read,n_search_total,&flag_go);
        if(flag_go){
           SID_log("Processing matches for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read);                 

           // Read base group
           sprintf(filename_cat1,"%03d",i_read);
           init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
           //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
           read_groups(filename_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,&plist1,filename_cat1,-1,-1);
           //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

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

                 // If this snapshot combination is missing at least one file, proceed.
                 flag_read=TRUE;
                 if(flag_go){
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
                          //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
                          read_groups(filename_root_in,j_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,
                                      &plist2,filename_cat2,PHK_min_local,PHK_max_local);
                          //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
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
                          //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
                          match_halos(plist1_order,i_read_order,NULL,0,plist2_order,j_read_order,NULL,0,"match",flag_match_subgroups|MATCH_STORE_SCORE);
                          //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);

                          // Writing results
                          //SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
                          write_match_results(filename_out_dir,filename_out_name,i_read_order,j_read_order,plist1_order,plist2_order,k_match);
                          //SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
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
        } // If this snapshot needs to be processed
     } // Loop over base snapshots
     SID_log("Done.",SID_LOG_CLOSE);
  }

  // Read the stuff that needs to be returned back to the calling function
  SID_log("Reading header information...",SID_LOG_OPEN);

  // Count the number of snapshots we are going to use and
  //    initialize the arrays that are to be returned
  for(i_read=i_read_stop,(*n_files_return)=0;      
      i_read>=i_read_start;
      i_read-=i_read_step) (*n_files_return)++;  
  (*n_subgroups_return)=(int *)SID_malloc(sizeof(int)*(*n_files_return));
  (*n_groups_return)   =(int *)SID_malloc(sizeof(int)*(*n_files_return));

  if(SID.I_am_Master){
     FILE *fp_read_header;
     for(k_match=0;k_match<2;k_match++){
        switch(k_match){
           case 0:
           sprintf(group_text_prefix,"sub");
           flag_compute_header=flag_compute_header_subgroups;
           break;
           case 1:
           sprintf(group_text_prefix,"");
           flag_compute_header=flag_compute_header_groups;
           break;
        }
        for(i_read=i_read_stop,j_read=0;i_read>=i_read_start;i_read-=i_read_step,j_read++){
           sprintf(filename_cat1,"%03d",i_read);
           sprintf(filename_out,"%s/%s.%sgroup_matches_header",filename_out_dir,filename_out_name,group_text_prefix);

           // Open file and skip header           
           if(i_read==i_read_stop){
              if((fp_read_header=fopen(filename_out,"r"))==NULL)
                 SID_trap_error("Could not open file {%s} when reading header information.",ERROR_IO_OPEN,filename_out);
              fseek(fp_read_header,4*sizeof(int),SEEK_SET);
           }

           fseek(fp_read_header,1*sizeof(int),SEEK_CUR);
           fread(&n_groups_1,sizeof(int),1,fp_read_header);
           fseek(fp_read_header,n_groups_1*sizeof(int),SEEK_CUR);
           if(k_match==1)
              fseek(fp_read_header,n_groups_1*sizeof(int),SEEK_CUR);
           switch(k_match){
              case 0:
              (*n_subgroups_return)[j_read]=n_groups_1;
              break;
              case 1:
              (*n_groups_return)[j_read]   =n_groups_1;
              break;
           }
        }
        fclose(fp_read_header);
     }
  }
  SID_Bcast((*n_subgroups_return),sizeof(int)*(*n_files_return),MASTER_RANK,SID.COMM_WORLD);
  SID_Bcast((*n_groups_return),   sizeof(int)*(*n_files_return),MASTER_RANK,SID.COMM_WORLD);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);
}

