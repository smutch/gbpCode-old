#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void compute_trees_matches(char   *filename_halo_root_in,
                           char   *filename_root_out,
                           int     i_read_start,
                           int     i_read_stop,
                           int     i_read_step,
                           int    *n_files_return,
                           int   **n_subgroups_return,
                           int   **n_groups_return,
                           int     n_search){
  char        filename_out[256];
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
  char        filename_groups[256];
  int         n_groups_1;
  int         n_groups_2;
  int         i_rank;
  int        *n_particles;
  int        *n_sub_group;
  int        *match_id   =NULL;
  float      *match_score=NULL;
  char        cat_name_1[20];
  char        cat_name_2[20];
  size_t     *match_index=NULL;
  size_t      offset;
  plist_info  plist1;
  plist_info  plist2;
  int        *n_return;

  SID_log("Constructing merger tree matches...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Initialize a few things
  init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  init_plist(&plist2,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
  n_search_total=n_search*i_read_step;

  // Initialize the n_halos arrays that are to be returned
  for(i_read=i_read_stop,(*n_files_return)=0;      
      i_read>=i_read_start;
      i_read-=i_read_step) (*n_files_return)++;  
  (*n_subgroups_return)=(int *)SID_malloc(sizeof(int)*(*n_files_return));
  (*n_groups_return)   =(int *)SID_malloc(sizeof(int)*(*n_files_return));

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
     sprintf(filename_out,"%s.%sgroup_matches",filename_root_out,group_text_prefix);

     flag_go=TRUE;
     if((fp_test=fopen(filename_out,"r"))!=NULL){
        fclose(fp_test);
        SID_log("Checking if existing %sgroup matching file is adequate...",SID_LOG_OPEN,group_text_prefix);
        SID_fopen(filename_out,"r",&fp_in);
        SID_fread(&i_read_start_file,sizeof(int),1,&fp_in);
        SID_fread(&i_read_stop_file, sizeof(int),1,&fp_in);
        SID_fread(&n_search_file,    sizeof(int),1,&fp_in);
        SID_fread(&n_files,          sizeof(int),1,&fp_in);
        SID_fclose(&fp_in);
        if(i_read_stop_file>=i_read_stop && i_read_start_file<=i_read_start && n_search_file>=n_search_total)
           flag_go=FALSE;
        else{
           flag_go=TRUE;
           SID_log("Matches will need to be recomputed.",SID_LOG_COMMENT);
        }
        SID_log("Done.",SID_LOG_CLOSE);
     }
     if(flag_go){
        SID_log("Constructing %sgroup merger tree matches...",SID_LOG_OPEN,group_text_prefix);
        SID_log("Writing header...",SID_LOG_OPEN);
        SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
        SID_fopen(filename_out,"w",&fp_out);
        offset=0;
        SID_fwrite(&i_read_start,  sizeof(int),1,&fp_out);offset+=sizeof(int);
        SID_fwrite(&i_read_stop,   sizeof(int),1,&fp_out);offset+=sizeof(int);
        SID_fwrite(&n_search_total,sizeof(int),1,&fp_out);offset+=sizeof(int);
        // Count the number of matches and files that we will store
        for(i_read=i_read_stop,n_matches=0,n_files=0;i_read>=i_read_start;i_read--){
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 n_matches++;
              }
           }
           n_files++;
        }
        // Write halo counts for each file
        SID_fwrite(&n_files,sizeof(int),1,&fp_out);offset+=sizeof(int);
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           sprintf(filename_cat1,  "%03d",i_read);
           sprintf(filename_groups,"%s_%03d",filename_halo_root_in,i_read,group_text_prefix);
           read_groups(filename_halo_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist1,filename_cat1);
           n_groups_1 =((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_%s", group_text_prefix,filename_cat1))[0];
           SID_fwrite(&i_read,    sizeof(int),1,         &fp_out);offset+=sizeof(int);
           SID_fwrite(&n_groups_1,sizeof(int),1,         &fp_out);offset+=sizeof(int);
           if(n_groups_1>0){
              n_particles= (int    *)ADaPS_fetch(plist1.data,"n_particles_%sgroup_%s",group_text_prefix,filename_cat1);
              n_sub_group= (int    *)ADaPS_fetch(plist1.data,"n_subgroups_group_%s",filename_cat1);
              SID_fwrite(n_particles,sizeof(int),n_groups_1,&fp_out);offset+=sizeof(int)*n_groups_1;
              if(flag_match_subgroups==MATCH_GROUPS){
                 SID_fwrite(n_sub_group,sizeof(int),n_groups_1,&fp_out);offset+=sizeof(int)*n_groups_1;
              }
           }
           free_plist(&plist1);
        }
        SID_fwrite(&n_matches,sizeof(int),1,&fp_out);offset+=sizeof(int);
        // Finish computing the offset to the start of the matching data
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 offset+=2*sizeof(int)+sizeof(size_t);
              }
           }
        }
        // Write matching combinations and the offsets to each's matching data
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 sprintf(filename_cat1,  "%03d",i_read);
                 sprintf(filename_groups,"%s_%03d",filename_halo_root_in,i_read,group_text_prefix);
                 read_groups(filename_halo_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_NOIDS,&plist1,filename_cat1);
                 n_groups_1 =((int *)ADaPS_fetch(plist1.data,"n_%sgroups_%s", group_text_prefix,filename_cat1))[0];
                 SID_fwrite(&i_read,sizeof(int),   1,&fp_out);
                 SID_fwrite(&j_read,sizeof(int),   1,&fp_out);
                 SID_fwrite(&offset,sizeof(size_t),1,&fp_out);
                 offset+=(4*sizeof(int)+n_groups_1*(sizeof(int)+sizeof(size_t)+sizeof(float)));
                 free_plist(&plist1);
              }
           }
        }
        SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
        SID_log("Done.",SID_LOG_CLOSE);
        // Compute each matching combaination and write the results to the file
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
             SID_log("Processing %sgroup matches for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix,i_read);                 
             SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

             // Read base group
             sprintf(filename_cat1,  "%03d",i_read);
             sprintf(filename_groups,"%s_%03d",filename_halo_root_in,i_read,group_text_prefix);
             read_groups(filename_halo_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist1,filename_cat1);
             n_groups_1=((int *)ADaPS_fetch(plist1.data,"n_%sgroups_%s",group_text_prefix,filename_cat1))[0];
                      
             for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
                if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                   SID_log("Processing %sgroup matches to snap %03d...",SID_LOG_OPEN,group_text_prefix,j_read);                 

                   // Read catalog to match to
                   sprintf(filename_cat2,  "%03d",j_read);
                   sprintf(filename_groups,"%s_%03d",filename_halo_root_in,j_read,group_text_prefix);
                   read_groups(filename_halo_root_in,j_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist2,filename_cat2);

                   // Perform matching
                   match_halos(&plist1,i_read,NULL,0,&plist2,j_read,NULL,0,"match",flag_match_subgroups|MATCH_STORE_SCORE);

                   // Write results to output
                   for(i_rank=0;i_rank<SID.n_proc;i_rank++){
                      if(SID.My_rank==i_rank){
                        // Write a few simple stats at start
                         if(i_rank==0){
                            match_id   =(int   *)ADaPS_fetch(plist1.data,"match_match");
                            match_score=(float *)ADaPS_fetch(plist1.data,"match_score_match");
                            n_groups_2 =((int *)ADaPS_fetch(plist2.data,"n_%sgroups_%s",group_text_prefix,filename_cat2))[0];
                            merge_sort((void *)match_id,(size_t)n_groups_1,&match_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
                            SID_fwrite(&i_read,    sizeof(int),   1,         &fp_out);
                            SID_fwrite(&j_read,    sizeof(int),   1,         &fp_out);
                            SID_fwrite(&n_groups_1,sizeof(int),   1,         &fp_out);
                            SID_fwrite(&n_groups_2,sizeof(int),   1,         &fp_out);
                            SID_fwrite(match_id,   sizeof(int),   n_groups_1,&fp_out);
                            SID_fwrite(match_index,sizeof(size_t),n_groups_1,&fp_out);
                            SID_fwrite(match_score,sizeof(float), n_groups_1,&fp_out);
                            SID_free(SID_FARG match_index);
                         }
                      }
                   }
                   free_plist(&plist2);
                   SID_log("Done.",SID_LOG_CLOSE);
                }
             }
             free_plist(&plist1);
             SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
             SID_log("Done.",SID_LOG_CLOSE);
        }
        SID_fclose(&fp_out);
        SID_log("Done.",SID_LOG_CLOSE);
     }

     // Construct the n_halos arrays that are to be returned
     SID_fopen(filename_out,"r",&fp_in);
     SID_fread(&i_read_start_file, sizeof(int),1,&fp_in);
     SID_fread(&i_read_stop_file,sizeof(int),1,&fp_in);
     SID_fread(&n_search,         sizeof(int),1,&fp_in);
     SID_fread(&n_files,          sizeof(int),1,&fp_in);
     for(i_read=i_read_stop_file,j_read=i_read,k_read=0;i_read>=i_read_start_file;i_read--){
        SID_fread(&l_read,    sizeof(int),1,         &fp_in);
        SID_fread(&n_groups_1,sizeof(int),1,         &fp_in);
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
     if(k_read!=(*n_files_return)) SID_trap_error("%sgroup size arrays not properly populated (i.e. %d!=%d).",ERROR_LOGIC,group_text_prefix,k_read,(*n_files_return));
  }
  free_plist(&plist2);
  free_plist(&plist1);
  SID_log("Done.",SID_LOG_CLOSE);
  
}
