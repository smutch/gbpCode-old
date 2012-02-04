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

  SID_log("Constructing merger tree matches...",SID_LOG_OPEN|SID_LOG_TIMER);

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
        SID_fread_all(&i_read_start_file,sizeof(int),1,&fp_in);
        SID_fread_all(&i_read_stop_file, sizeof(int),1,&fp_in);
        SID_fread_all(&n_search_file,    sizeof(int),1,&fp_in);
        SID_fread_all(&n_files,          sizeof(int),1,&fp_in);
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
        // Write the header
        SID_log("Writing header...",SID_LOG_OPEN|SID_LOG_TIMER);
        SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
        // ... write snapshot ranges and search range ...
        offset=0;
        if(SID.I_am_Master){
           SID_fopen(filename_out,"w",&fp_out);
           SID_fwrite(&i_read_start,  sizeof(int),1,&fp_out);offset+=sizeof(int);
           SID_fwrite(&i_read_stop,   sizeof(int),1,&fp_out);offset+=sizeof(int);
           SID_fwrite(&n_search_total,sizeof(int),1,&fp_out);offset+=sizeof(int);
        }

        // Count the number of matches and files that we will store
        for(i_read=i_read_stop,n_matches=0,n_files=0;i_read>=i_read_start;i_read--){
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 n_matches++;
              }
           }
           n_files++;
        }
        // ... write the number of files involved ...
        if(SID.I_am_Master){
           SID_fwrite(&n_files,sizeof(int),1,&fp_out);offset+=sizeof(int);
        }

        // Write halo counts for each file in the given range
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           sprintf(filename_cat1,"%03d",i_read);
           init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
           read_groups(filename_halo_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_NOIDS,&plist1,filename_cat1);
           n_groups_1      =((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_all_%s",group_text_prefix,filename_cat1))[0];
           n_groups_1_local=((int    *)ADaPS_fetch(plist1.data,"n_%sgroups_%s",    group_text_prefix,filename_cat1))[0];
           if(SID.I_am_Master){
              SID_fwrite(&i_read,    sizeof(int),1,&fp_out);offset+=sizeof(int);
              SID_fwrite(&n_groups_1,sizeof(int),1,&fp_out);offset+=sizeof(int);
           }
           if(n_groups_1>0){
              n_particles=(int *)ADaPS_fetch(plist1.data,"n_particles_%sgroup_%s",group_text_prefix,filename_cat1);
              SID_fwrite_ordered(n_particles,sizeof(int),(size_t)n_groups_1_local,&fp_out);offset+=n_groups_1*sizeof(int);
              if(flag_match_subgroups==MATCH_GROUPS){
                 n_sub_group=(int *)ADaPS_fetch(plist1.data,"n_subgroups_group_%s",filename_cat1);
                 SID_fwrite_ordered(n_sub_group,sizeof(int),(size_t)n_groups_1_local,&fp_out);offset+=n_groups_1*sizeof(int);
              }
           }
           free_plist(&plist1);
        }
        if(SID.I_am_Master){
           SID_fwrite(&n_matches,sizeof(int),1,&fp_out);offset+=sizeof(int);
        }

        // Finish computing the offset to the start of the matching data.  These are due to
        //   the writes in the nested loop following this one.
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 offset+=2*sizeof(int)+sizeof(size_t);
              }
           }
        }
        // Write matching combinations and the offsets to each's matching data
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           sprintf(filename_cat1,  "%03d",i_read);
           init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
           read_groups(filename_halo_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_NOIDS,&plist1,filename_cat1);
           n_groups_1      =((int *)ADaPS_fetch(plist1.data,"n_%sgroups_all_%s",group_text_prefix,filename_cat1))[0];
           n_groups_1_local=((int *)ADaPS_fetch(plist1.data,"n_%sgroups_%s",    group_text_prefix,filename_cat1))[0];
           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 if(SID.I_am_Master){
                    SID_fwrite(&i_read,sizeof(int),   1,&fp_out);
                    SID_fwrite(&j_read,sizeof(int),   1,&fp_out);
                    SID_fwrite(&offset,sizeof(size_t),1,&fp_out);
                 }
                 offset+=(4*sizeof(int)+n_groups_1*(sizeof(int)+sizeof(size_t)+sizeof(float)));
              }
           }
           free_plist(&plist1);
        }

        SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
        SID_log("Done.",SID_LOG_CLOSE);

        // Compute each matching combaination and write the results to the file
        for(i_read=i_read_stop;i_read>=i_read_start;i_read--){
           SID_log("Processing %sgroup matches for snapshot #%d...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix,i_read);                 
           SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);

           // Read base group
           int *file_index_1;
           int  PHK_min_local;
           int  PHK_max_local;
           sprintf(filename_cat1,  "%03d",i_read);
           init_plist(&plist1,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
           read_groups(filename_halo_root_in,i_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,&plist1,filename_cat1,-1,-1);
           PHK_min_local   =((int *)ADaPS_fetch(plist1.data,"PHK_min_local_%s",      filename_cat1))[0];
           PHK_max_local   =((int *)ADaPS_fetch(plist1.data,"PHK_max_local_%s",      filename_cat1))[0];
           n_groups_1      =((int *)ADaPS_fetch(plist1.data,"n_%sgroups_all_%s",     group_text_prefix,filename_cat1))[0];
           n_groups_1_local=((int *)ADaPS_fetch(plist1.data,"n_%sgroups_%s",         group_text_prefix,filename_cat1))[0];
           file_index_1    = (int *)ADaPS_fetch(plist1.data,"file_index_%sgroups_%s",group_text_prefix,filename_cat1);
// Test read_groups
/*
FILE   *fp_test;
char    filename_test[256];
int    *group_size;
int    *group_offset;
size_t *ids;
int     buffer_0[10000];
sprintf(filename_test,"test_%sgroup_size_%s.dat",group_text_prefix,filename_cat1);
group_size  =(int    *)ADaPS_fetch(plist1.data,"n_particles_%sgroup_%s",group_text_prefix,filename_cat1);
group_offset=(int    *)ADaPS_fetch(plist1.data,"particle_offset_%sgroup_%s",group_text_prefix,filename_cat1);
ids         =(size_t *)ADaPS_fetch(plist1.data,"particle_ids_%s",filename_cat1);
fp_test=fopen(filename_test,"w");
for(i_group=0,buffered_count_local=0;i_group<n_groups_1;i_group+=n_buffer){
   // Decide this buffer iteration's size
   n_buffer=MIN(n_buffer_max,n_groups_1-i_group);
   // Set the buffer to a default value smaller than the smallest possible data size
   for(i_buffer=0;i_buffer<n_buffer;i_buffer++)
      buffer_0[i_buffer]=-2; // Min value of match_id is -1
   // Determine if any of the local data is being used for this buffer
   for(j_group=0;j_group<n_groups_1_local;j_group++){
      index_test=file_index_1[j_group]-i_group;
      // ... if so, set the appropriate buffer value
      if(index_test>=0 && index_test<n_buffer){
        buffer_0[index_test]=group_size[j_group];
        buffered_count_local++;
      }
   }
   // Doing a global max on the buffer yields the needed buffer on all ranks
   SID_Allreduce(SID_IN_PLACE,buffer_0,n_buffer,SID_INT,SID_MAX,SID.COMM_WORLD);
   // Write the buffer
   if(SID.I_am_Master){
     for(i_buffer=0;i_buffer<n_buffer;i_buffer++) fprintf(fp_test,"%6d %6d\n",i_group+i_buffer,buffer_0[i_buffer]);
   }
}
fclose(fp_test);
sprintf(filename_test,"test_%sgroup_IDs_%s.dat",group_text_prefix,filename_cat1);
fp_test=fopen(filename_test,"w");
size_t buffer_1[10000];
size_t buffer_2[10000];
for(i_group=0,buffered_count_local=0;i_group<n_groups_1;i_group+=n_buffer){
   // Decide this buffer iteration's size
   n_buffer=MIN(n_buffer_max,n_groups_1-i_group);
   // Set the buffer to a default value smaller than the smallest possible data size
   for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
      buffer_1[i_buffer]=0;
      buffer_2[i_buffer]=0;
   }
   // Determine if any of the local data is being used for this buffer
   for(j_group=0;j_group<n_groups_1_local;j_group++){
      index_test=file_index_1[j_group]-i_group;
      // ... if so, set the appropriate buffer value
      if(index_test>=0 && index_test<n_buffer){
        if(group_size[j_group]>0){
          buffer_1[index_test]=ids[group_offset[j_group]];
          buffer_2[index_test]=ids[group_offset[j_group]+group_size[j_group]-1];
        }
        buffered_count_local++;
      }
   }
   // Doing a global max on the buffer yields the needed buffer on all ranks
   SID_Allreduce(SID_IN_PLACE,buffer_1,n_buffer,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
   SID_Allreduce(SID_IN_PLACE,buffer_2,n_buffer,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
   // Write the buffer
   if(SID.I_am_Master){
     for(i_buffer=0;i_buffer<n_buffer;i_buffer++) fprintf(fp_test,"%6d %6lld %6lld\n",i_group+i_buffer,buffer_1[i_buffer],buffer_2[i_buffer]);
   }
}
fclose(fp_test);
SID_exit(ERROR_NONE);
*/

           for(j_read=i_read+n_search_total;j_read>=i_read-n_search_total;j_read--){
              if(j_read<=i_read_stop && j_read>=i_read_start && j_read!=i_read){
                 SID_log("Processing %sgroup matches to snap %03d...",SID_LOG_OPEN,group_text_prefix,j_read);                 

                 // Read catalog to match to
                 sprintf(filename_cat2,  "%03d",j_read);
                 init_plist(&plist2,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
                 read_groups(filename_halo_root_in,j_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES|READ_GROUPS_PEANOHILBERT,&plist2,filename_cat2,PHK_min_local,PHK_max_local);
                 //read_groups(filename_halo_root_in,j_read,READ_GROUPS_ALL|READ_GROUPS_NOPROPERTIES,&plist2,filename_cat2);

                 // Perform matching
                 match_halos(&plist1,i_read,NULL,0,&plist2,j_read,NULL,0,"match",flag_match_subgroups|MATCH_STORE_SCORE);

                 // Write results to output ...
                 match_id        = (int   *)ADaPS_fetch(plist1.data,"match_match");
                 match_score     = (float *)ADaPS_fetch(plist1.data,"match_score_match");
                 n_groups_2      =((int   *)ADaPS_fetch(plist2.data,"n_%sgroups_all_%s",group_text_prefix,filename_cat2))[0];
                 n_groups_2_local=((int   *)ADaPS_fetch(plist2.data,"n_%sgroups_%s",    group_text_prefix,filename_cat2))[0];
                 sort(match_id,   (size_t)n_groups_1_local,&match_index,SID_INT,   SORT_GLOBAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
                 sort(match_index,(size_t)n_groups_1_local,&match_rank, SID_SIZE_T,SORT_GLOBAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

                 // ... write a few simple stats at start ...
                 if(SID.I_am_Master){
                    SID_fwrite(&i_read,    sizeof(int),1,&fp_out);
                    SID_fwrite(&j_read,    sizeof(int),1,&fp_out);
                    SID_fwrite(&n_groups_1,sizeof(int),1,&fp_out);
                    SID_fwrite(&n_groups_2,sizeof(int),1,&fp_out);
                 }
                 SID_Barrier(SID.COMM_WORLD);

                 // Write matching results.  We need to write back to the file in the
                 //   order that it was read from the halo catalogs, not in the PH order
                 //   that it is stored in RAM.  This requires some buffering.
                 buffer       =SID_malloc(n_buffer_max*sizeof(size_t));
                 buffer_int   =(int    *)buffer;
                 buffer_size_t=(size_t *)buffer;
                 buffer_float =(float  *)buffer;

                 // Write match_ids ...
                 //    ... loop over all the groups in buffer-sized batches
FILE *fp_test;
char  filename_test[256];
size_t buffer_test[1000];
sprintf(filename_test,"test_%sgroup_match_IDs_%s_%s.dat",group_text_prefix,filename_cat1,filename_cat2);
if(SID.I_am_Master) fp_test=fopen(filename_test,"w");
                 for(i_group=0,buffered_count_local=0;i_group<n_groups_1;i_group+=n_buffer){
                    // Decide this buffer iteration's size
                    n_buffer=MIN(n_buffer_max,n_groups_1-i_group);
                    // Set the buffer to a default value smaller than the smallest possible data size
                    for(i_buffer=0;i_buffer<n_buffer;i_buffer++)
                       buffer_int[i_buffer]=-2; // Min value of match_id is -1
for(i_buffer=0;i_buffer<n_buffer;i_buffer++)
   buffer_test[i_buffer]=0;
                    // Determine if any of the local data is being used for this buffer
                    for(j_group=0;j_group<n_groups_1_local;j_group++){
                       index_test=file_index_1[j_group]-i_group;
                       // ... if so, set the appropriate buffer value
                       if(index_test>=0 && index_test<n_buffer){
                         buffer_int[index_test]=match_id[j_group];
buffer_test[index_test]=match_rank[j_group];
                         buffered_count_local++;
                       }
                    }
                    // Doing a global max on the buffer yields the needed buffer on all ranks
                    SID_Allreduce(SID_IN_PLACE,buffer_int,n_buffer,SID_INT,SID_MAX,SID.COMM_WORLD);
SID_Allreduce(SID_IN_PLACE,buffer_test,n_buffer,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
                    // Sanity check
                    for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
                      if(buffer_int[i_buffer]<-1)
                        SID_trap_error("Illegal match_id result (%d) for group No. %d.",ERROR_LOGIC,buffer_int[i_buffer],i_group+i_buffer);
                    }
                    // Write the buffer
                    if(SID.I_am_Master){
                      SID_fwrite(buffer_int,sizeof(int),(size_t)n_buffer,&fp_out);
for(i_buffer=0;i_buffer<n_buffer;i_buffer++) fprintf(fp_test,"%7d %7d %7lld\n",i_group+i_buffer,buffer_int[i_buffer],buffer_test[i_buffer]);
                    }
                 }
if(SID.I_am_Master) fclose(fp_test);
                 // Sanity check
                 calc_sum_global(&buffered_count_local,&buffered_count,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
                 if(buffered_count!=n_groups_1)
                   SID_trap_error("Buffer counts don't make sense (ie %d!=%d) after writing match IDs.",ERROR_LOGIC,buffered_count,n_groups_1);

                 // Write match sort indices ...
                 //    ... loop over all the groups in buffer-sized batches
sprintf(filename_test,"test_%sgroup_match_index_%s_%s.dat",group_text_prefix,filename_cat1,filename_cat2);
if(SID.I_am_Master) fp_test=fopen(filename_test,"w");
                 for(i_group=0,buffered_count_local=0;i_group<n_groups_1;i_group+=n_buffer){
                    // Decide this buffer iteration's size
                    n_buffer=MIN(n_buffer_max,n_groups_1-i_group);
                    // Set the buffer to a default value smaller than the smallest possible data size
                    for(i_buffer=0;i_buffer<n_buffer;i_buffer++)
                       buffer_size_t[i_buffer]=-1; // Min value of match_rank is 0
                    // Determine if any of the local data is being used for this buffer
                    for(j_group=0;j_group<n_groups_1_local;j_group++){
                       index_test=match_rank[j_group]-i_group;
                       // ... if so, set the appropriate buffer value
                       if(index_test>=0 && index_test<n_buffer){
                         buffer_size_t[index_test]=file_index_1[j_group];
                         buffered_count_local++;
                       }
                    }
                    // Doing a global max on the buffer yields the needed buffer on all ranks
                    SID_Allreduce(SID_IN_PLACE,buffer_size_t,n_buffer,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
                    // Sanity check
                    for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
                      if(buffer_size_t[i_buffer]<0)
                        SID_trap_error("Illegal match_rank result (%lld) for group No. %d.",ERROR_LOGIC,buffer_size_t[i_buffer],i_group+i_buffer);
                    }
                    // Write the buffer
                    if(SID.I_am_Master){
                      SID_fwrite(buffer,sizeof(size_t),(size_t)n_buffer,&fp_out);
for(i_buffer=0;i_buffer<n_buffer;i_buffer++) fprintf(fp_test,"%7d %7lld\n",i_group+i_buffer,buffer_size_t[i_buffer]);
                    }
                 }
if(SID.I_am_Master) fclose(fp_test);
                 // Sanity check
                 calc_sum_global(&buffered_count_local,&buffered_count,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
                 if(buffered_count!=n_groups_1)
                   SID_trap_error("Buffer counts don't make sense (ie %d!=%d) after writing match indices.",ERROR_LOGIC,buffered_count,n_groups_1);

                 // Write match_score ...
                 //    ... loop over all the groups in buffer-sized batches
                 for(i_group=0,buffered_count_local=0;i_group<n_groups_1;i_group+=n_buffer){
                    // Decide this buffer iteration's size
                    n_buffer=MIN(n_buffer_max,n_groups_1-i_group);
                    // Set the buffer to a default value smaller than the smallest possible data size
                    for(i_buffer=0;i_buffer<n_buffer;i_buffer++)
                       buffer_float[i_buffer]=-1.; // Min value of match_score is 0.
                    // Determine if any of the local data is being used for this buffer
                    for(j_group=0;j_group<n_groups_1_local;j_group++){
                       index_test=file_index_1[j_group]-i_group;
                       // ... if so, set the appropriate buffer value
                       if(index_test>=0 && index_test<n_buffer){
                         buffer_float[index_test]=match_score[j_group];
                         buffered_count_local++;
                       }
                    }
                    // Doing a global max on the buffer yields the needed buffer on all ranks
                    SID_Allreduce(SID_IN_PLACE,buffer_float,n_buffer,SID_FLOAT,SID_MAX,SID.COMM_WORLD);
                    // Sanity check
                    for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
                      if(buffer_float[i_buffer]<-0.)
                        SID_trap_error("Illegal match_score result (%f) for group No. %d.",ERROR_LOGIC,buffer_float[i_buffer],i_group+i_buffer);
                    }
                    // Write the buffer
                    if(SID.I_am_Master)
                      SID_fwrite(buffer,sizeof(float),(size_t)n_buffer,&fp_out);
                 }
                 // Sanity check
                 calc_sum_global(&buffered_count_local,&buffered_count,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
                 if(buffered_count!=n_groups_1)
                    SID_trap_error("Buffer counts don't make sense (ie %d!=%d) after writing match scores.",ERROR_LOGIC,buffered_count,n_groups_1);

                 // Clean-up
                 SID_free(SID_FARG buffer);
                 SID_free(SID_FARG match_rank);
                 free_plist(&plist2);
                 SID_log("Done.",SID_LOG_CLOSE);
              }
           }
           free_plist(&plist1);
           SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
           SID_log("Done.",SID_LOG_CLOSE);
        } // loop over snapshots
        if(SID.I_am_Master)
           SID_fclose(&fp_out);
        SID_log("Done.",SID_LOG_CLOSE);
     } // if(flag_go)

     // Construct the n_halos arrays that are to be returned
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
//     if(k_read!=(*n_files_return)) SID_trap_error("%sgroup size arrays not properly populated (i.e. %d!=%d).",ERROR_LOGIC,group_text_prefix,k_read,(*n_files_return));
  } // loop over k_match
  SID_log("Done.",SID_LOG_CLOSE);
}

