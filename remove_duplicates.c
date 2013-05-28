#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpHalos.h>

#define GADGET_BUFFER_SIZE_LOCAL  128*SIZE_OF_MEGABYTE

void read_duplicates_file_local(const char  *filename_in,const char *filename_number,
                                plist_info  *plist,
                                int         *size_group,
                                size_t      *offset_group,
                                int         *size_subgroup,
                                size_t      *offset_subgroup,
                                int         *n_subgroups_group,
                                int          n_groups_local,
                                int          n_subgroups_local,
                                int          group_index_local_offset,
                                int          subgroup_index_local_offset,
                                int         *flag_add_duplicates,
                                int         *n_duplicates_added,
                                int         *n_duplicates_added_local,
                                size_t     **id_duplicate_added,
                                int        **group_duplicate_added_to,
                                int        **subgroup_duplicate_added_to,
                                int        **index_duplicate_added_to);
void read_duplicates_file_local(const char  *filename_in,const char *filename_number,
                                plist_info  *plist,
                                int         *size_group,
                                size_t      *offset_group,
                                int         *size_subgroup,
                                size_t      *offset_subgroup,
                                int         *n_subgroups_group,
                                int          n_groups_local,
                                int          n_subgroups_local,
                                int          group_index_local_offset,
                                int          subgroup_index_local_offset,
                                int         *flag_add_duplicates,
                                int         *n_duplicates_added,
                                int         *n_duplicates_added_local,
                                size_t     **id_duplicate_added,
                                int        **group_duplicate_added_to,
                                int        **subgroup_duplicate_added_to,
                                int        **index_duplicate_added_to){
   FILE *fp_in;
   if((fp_in=fopen(filename_in,"r"))==NULL)
      (*flag_add_duplicates)=FALSE;
   else{
      SID_log("DUPLICATE TEST FILE FOUND.",SID_LOG_COMMENT);
      size_t *id_list_local;
      size_t  n_particles_all;
      size_t  n_particles_local;
      n_particles_all  =((size_t *)ADaPS_fetch(plist->data,"n_particles_all_%s",filename_number))[0];
      n_particles_local=((size_t *)ADaPS_fetch(plist->data,"n_particles_%s",    filename_number))[0];
      id_list_local    = (size_t *)ADaPS_fetch(plist->data,"particle_ids_%s",   filename_number);
      size_t line_length=0;
      char  *line       =NULL;
      int group_duplicate_added_from;
      int index_duplicate_added_from;
      (*n_duplicates_added)         =count_lines_data(fp_in);
      (*id_duplicate_added)         =(size_t *)SID_calloc(sizeof(size_t)*(*n_duplicates_added));
      (*group_duplicate_added_to)   =(int    *)SID_malloc(sizeof(int)   *(*n_duplicates_added));
      (*subgroup_duplicate_added_to)=(int    *)SID_calloc(sizeof(int)   *(*n_duplicates_added));
      (*index_duplicate_added_to)   =(int    *)SID_malloc(sizeof(int)   *(*n_duplicates_added));
      int i_duplicate;
      int index_test;
      (*n_duplicates_added_local)=0;
      for(i_duplicate=0;i_duplicate<(*n_duplicates_added);i_duplicate++){
         grab_next_line_data(fp_in,&line,&line_length);
         grab_int(line,1,&group_duplicate_added_from);
         grab_int(line,2,&index_duplicate_added_from);
         grab_int(line,3,&((*group_duplicate_added_to)[i_duplicate]));
         grab_int(line,4,&((*index_duplicate_added_to)[i_duplicate]));
         // Set the particle IDs of the injected duplicates
         index_test=group_duplicate_added_from-group_index_local_offset;
         if(index_test>=0 && index_test<n_groups_local){
            (*n_duplicates_added_local)++;
            if(index_duplicate_added_from<0 || index_duplicate_added_from>=size_group[index_test])
               SID_trap_error("Duplicate No. %d's origin index is out of the specified group's index range.",ERROR_LOGIC,i_duplicate);
            index_duplicate_added_from+=offset_group[index_test];
            (*id_duplicate_added)[i_duplicate]=id_list_local[index_duplicate_added_from];
         }
         // Set the subgroup destination index
         index_test=(*group_duplicate_added_to)[i_duplicate]-group_index_local_offset;
         if(index_test>=0 && index_test<n_groups_local){
            int i_group;
            int i_subgroup;
            if((*index_duplicate_added_to)[i_duplicate]<0 || (*index_duplicate_added_to)[i_duplicate]>size_group[index_test])
               SID_trap_error("Duplicate No. %d's destination index is out of the specified group's index range.",ERROR_LOGIC,i_duplicate);
            (*index_duplicate_added_to)[i_duplicate]+=offset_group[index_test]; // Temporarily place things in the local ID list frame
            for(i_group=0,i_subgroup=0;i_group<index_test;i_group++)
               i_subgroup+=n_subgroups_group[i_group];
            while((offset_subgroup[i_subgroup]+size_subgroup[i_subgroup])<(*index_duplicate_added_to)[i_duplicate] && i_subgroup<(n_subgroups_local-1))
               i_subgroup++;
            if((*index_duplicate_added_to)[i_duplicate]>=offset_subgroup[i_subgroup] &&
               (*index_duplicate_added_to)[i_duplicate]<(offset_subgroup[i_subgroup]+size_subgroup[i_subgroup]))
               (*subgroup_duplicate_added_to)[i_duplicate]=i_subgroup;
            else
               (*subgroup_duplicate_added_to)[i_duplicate]=-1;
            (*index_duplicate_added_to)[i_duplicate]-=offset_group[index_test]; // Place things back into the group ID frame to avoid 
                                                                                //    complications from other injected particles
         }
         // This duplicate will not be added locally.  Mark it as such.
         else{
            (*group_duplicate_added_to)[i_duplicate]   =-1;
            (*subgroup_duplicate_added_to)[i_duplicate]=-1;
            (*index_duplicate_added_to)[i_duplicate]   =n_particles_local+1;
         }
      }
      SID_Allreduce(SID_IN_PLACE,(*id_duplicate_added),(*n_duplicates_added),SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
      fclose(fp_in);
      SID_free(SID_FARG line);
      (*flag_add_duplicates)=TRUE;
   } 
}

void inject_duplicates_local(plist_info *plist,const char *filename_number,
                             int     *size_group,
                             size_t  *offset_group,
                             int     *size_subgroup,
                             size_t  *offset_subgroup,
                             int      n_groups_local,
                             int      n_subgroups_local,
                             int      group_index_local_offset,
                             int      subgroup_index_local_offset,
                             int      n_duplicates_added,
                             int      n_duplicates_added_local,
                             size_t  *id_duplicate_added,
                             int     *group_duplicate_added_to,
                             int     *subgroup_duplicate_added_to,
                             int     *index_duplicate_added_to);
void inject_duplicates_local(plist_info *plist,const char *filename_number,
                             int     *size_group,
                             size_t  *offset_group,
                             int     *size_subgroup,
                             size_t  *offset_subgroup,
                             int      n_groups_local,
                             int      n_subgroups_local,
                             int      group_index_local_offset,
                             int      subgroup_index_local_offset,
                             int      n_duplicates_added,
                             int      n_duplicates_added_local,
                             size_t  *id_duplicate_added,
                             int     *group_duplicate_added_to,
                             int     *subgroup_duplicate_added_to,
                             int     *index_duplicate_added_to){
   // Loop over the duplicates, adding the specified test duplicates
   SID_log("Injecting %d test duplicates into catalog...",SID_LOG_OPEN|SID_LOG_TIMER,n_duplicates_added);

   size_t *id_list_local_in;
   size_t  n_particles_all;
   size_t  n_particles_local;
   n_particles_all  =((size_t *)ADaPS_fetch(plist->data,"n_particles_all_%s",filename_number))[0];
   n_particles_local=((size_t *)ADaPS_fetch(plist->data,"n_particles_%s",    filename_number))[0];
   id_list_local_in = (size_t *)ADaPS_fetch(plist->data,"particle_ids_%s",   filename_number);

   // Increase the local id array size
   size_t  n_realloc;
   size_t  n_realloc_global;
   size_t *id_list_local;
   n_realloc       =n_particles_local+(size_t)n_duplicates_added_local;
   n_realloc_global=n_particles_all  +(size_t)n_duplicates_added;
   id_list_local   =(size_t *)SID_calloc(sizeof(size_t)*n_realloc);
   memcpy(id_list_local,id_list_local_in,sizeof(size_t)*n_particles_local);
   ADaPS_store(&(plist->data),(void *)(&n_realloc),       "n_particles_%s",    ADaPS_SCALAR_SIZE_T,filename_number);
   ADaPS_store(&(plist->data),(void *)(&n_realloc_global),"n_particles_all_%s",ADaPS_SCALAR_SIZE_T,filename_number);
   ADaPS_store(&(plist->data),(void *)(id_list_local),    "particle_ids_%s",   ADaPS_DEFAULT,filename_number);

   int i_duplicate;
   int j_duplicate;
   int i_group;
   for(i_duplicate=0;i_duplicate<n_duplicates_added;i_duplicate++){
      int i_group;
      int i_subgroup;
      int j_subgroup;

      // Check if this duplicate belongs to a local group
      int index_i;
      index_i=group_duplicate_added_to[i_duplicate]-group_index_local_offset;
      if(index_i>=0 && index_i<n_groups_local){
         // Adjust subgroup information (if the duplicate belongs to one)
         if(subgroup_duplicate_added_to[i_duplicate]>=0){
            int i_subgroup;
            int index_j;
            index_j=subgroup_duplicate_added_to[i_duplicate];
            if(index_j<0 || index_j>=n_subgroups_local)
               SID_trap_error("The subgroup that duplicate No. %d is supposed to be added to is on the wrong rank.",ERROR_SYNTAX,i_duplicate);
            for(i_subgroup=(index_j+1);i_subgroup<n_subgroups_local;i_subgroup++)
               offset_subgroup[i_subgroup]++;
            size_subgroup[index_j]++;
         }

         // Adjust group information
         for(i_group=(index_i+1);i_group<n_groups_local;i_group++)
            offset_group[i_group]++;
         size_group[index_i]++;

         // Insert the ID into the local ID list
         size_t i_particle;
         size_t index_k;
         index_k=index_duplicate_added_to[i_duplicate]+offset_group[index_i];
         for(i_particle=n_particles_local;i_particle>index_k && i_particle>0;i_particle--)
            id_list_local[i_particle]=id_list_local[i_particle-1];
         if(i_particle!=index_k)
            SID_trap_error("Could not find the index needed to inject a duplicate.",ERROR_LOGIC);
         id_list_local[i_particle]=id_duplicate_added[i_duplicate];
         n_particles_local++;
      }
   }
   SID_log("Done.",SID_LOG_CLOSE);
}

void write_duplicates_local(const char *filename_out_root,
                            const char *filename_number,
                            int         n_duplicates_local,
                            size_t     *id_duplicates_local,
                            int        *index_duplicates_group_local,
                            int        *index_duplicates_subgroup_local,
                            GBPREAL    *r2_duplicates_local,
                            char       *flag_keep_local,
                            int         group_index_local_offset,
                            int         subgroup_index_local_offset);
void write_duplicates_local(const char *filename_out_root,
                            const char *filename_number,
                            int         n_duplicates_local,
                            size_t     *id_duplicates_local,
                            int        *index_duplicates_group_local,
                            int        *index_duplicates_subgroup_local,
                            GBPREAL    *r2_duplicates_local,
                            char       *flag_keep_local,
                            int         group_index_local_offset,
                            int         subgroup_index_local_offset){
   int      i_rank;
   int      i_duplicate;
   int      n_buffer;
   int      n_buffer_alloc;
   size_t  *buffer_ids;
   int     *buffer_group;
   int     *buffer_subgroup;
   GBPREAL *buffer_r2;
   char    *buffer_flag;
   char     filename_out[MAX_FILENAME_LENGTH];

   sprintf(filename_out,"%s_%s_duplicates.txt",filename_out_root,filename_number);
   SID_log("Writing duplicates to {%s}...",SID_LOG_OPEN,filename_out);

   // Allocate exchange buffers
   n_buffer_alloc=n_duplicates_local;
   SID_Allreduce(SID_IN_PLACE,&n_buffer_alloc,1,SID_INT,SID_MAX,SID.COMM_WORLD);
   buffer_ids     =(size_t  *)SID_malloc(sizeof(size_t)* n_buffer_alloc);
   buffer_group   =(int     *)SID_malloc(sizeof(int)*    n_buffer_alloc);
   buffer_subgroup=(int     *)SID_malloc(sizeof(int)*    n_buffer_alloc);
   buffer_r2      =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_buffer_alloc);
   buffer_flag    =(char    *)SID_malloc(sizeof(char)*   n_buffer_alloc);

   // Write header
   FILE *fp_out;
   if(SID.I_am_Master){
      int i_column=1;
      fp_out=fopen(filename_out,"w");
      fprintf(fp_out,"# Duplicate list for catalog {%s_%s}\n",filename_out_root,filename_number);
      fprintf(fp_out,"# Column (%02d): Duplicate ID\n",                                 i_column++);
      fprintf(fp_out,"#        (%02d): Duplicate host group\n",                         i_column++);
      fprintf(fp_out,"#        (%02d): Duplicate host subgroup (-ve means undefined)\n",i_column++);
      fprintf(fp_out,"#        (%02d): Duplicate keep flag (1=kept, 0=removed)\n",      i_column++);
      fprintf(fp_out,"#        (%02d): Duplicate MBP separation\n",                     i_column++);
   }

   // Perform buffer exchanges and writes
   for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      n_buffer=n_duplicates_local;
      memcpy(buffer_ids,     id_duplicates_local,            n_duplicates_local*sizeof(size_t));
      memcpy(buffer_group,   index_duplicates_group_local,   n_duplicates_local*sizeof(int));
      memcpy(buffer_subgroup,index_duplicates_subgroup_local,n_duplicates_local*sizeof(int));
      memcpy(buffer_r2,      r2_duplicates_local,            n_duplicates_local*sizeof(GBPREAL));
      memcpy(buffer_flag,    flag_keep_local,                n_duplicates_local*sizeof(char));
      for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++){
         buffer_group[i_duplicate]   +=group_index_local_offset;
         buffer_subgroup[i_duplicate]+=subgroup_index_local_offset;
      }
      SID_Bcast(&n_buffer,sizeof(int),i_rank,SID.COMM_WORLD);
      if(n_buffer>0){
         SID_Bcast(buffer_ids,     n_buffer*sizeof(size_t), i_rank,SID.COMM_WORLD);
         SID_Bcast(buffer_group,   n_buffer*sizeof(int),    i_rank,SID.COMM_WORLD);
         SID_Bcast(buffer_subgroup,n_buffer*sizeof(int),    i_rank,SID.COMM_WORLD);
         SID_Bcast(buffer_r2,      n_buffer*sizeof(GBPREAL),i_rank,SID.COMM_WORLD);
         SID_Bcast(buffer_flag,    n_buffer*sizeof(char),   i_rank,SID.COMM_WORLD);
         if(SID.I_am_Master){
            int i_buffer;
            for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
               fprintf(fp_out,"%10lld %7d %7d %2d %le\n",
                              buffer_ids[i_buffer],
                              buffer_group[i_buffer],
                              buffer_subgroup[i_buffer],
                              buffer_flag[i_buffer],
                              sqrt(buffer_r2[i_buffer]));
            }
         }
      }
   }
   if(SID.I_am_Master)
     fclose(fp_out);
   SID_free(SID_FARG buffer_ids);
   SID_free(SID_FARG buffer_group);
   SID_free(SID_FARG buffer_subgroup);
   SID_free(SID_FARG buffer_r2);
   SID_free(SID_FARG buffer_flag);
   SID_log("Done.",SID_LOG_CLOSE);
}

void construct_duplicate_list_local(plist_info *plist,
                                    int         snapshot_number,
                                    size_t     *id_list_local,
                                    size_t      n_particles_local,
                                    size_t      n_particles_all_in_groups,
                                    size_t    **id_duplicates_local_return,
                                    size_t    **index_duplicates_local_return,
                                    size_t    **index_duplicates_local_index_return,
                                    int       **index_duplicates_group_local_return,
                                    int       **index_duplicates_subgroup_local_return,
                                    int        *n_duplicates_global_return,
                                    int        *n_duplicates_local_return);
void construct_duplicate_list_local(plist_info *plist,
                                    int         snapshot_number,
                                    size_t     *id_list_local,
                                    size_t      n_particles_local,
                                    size_t      n_particles_all_in_groups,
                                    size_t    **id_duplicates_local_return,
                                    size_t    **index_duplicates_local_return,
                                    size_t    **index_duplicates_local_index_return,
                                    int       **index_duplicates_group_local_return,
                                    int       **index_duplicates_subgroup_local_return,
                                    int        *n_duplicates_global_return,
                                    int        *n_duplicates_local_return){

  // Sort IDs globally
  size_t *id_list_local_rank =NULL;
  sort(id_list_local,n_particles_local,&id_list_local_rank,SID_SIZE_T,SORT_GLOBAL,SORT_COMPUTE_RANK, SORT_COMPUTE_NOT_INPLACE); 

  // Obtain a list of sorted IDs for each processor one buffer at a time
  size_t *buffer_id;
  size_t *buffer_index_local;
  char   *buffer_flag_local;
  int     n_buffer_max=1024*1024;
  int     n_buffer;
  size_t  id_last;
  size_t  index_last;
  char    flag_last;
  int     i_buffer;
  int     index_rank_test_local;  
  int     n_duplicates_unique_global=0;
  int     n_duplicates_global=0;
  int     n_duplicates_local =0;
  size_t *id_duplicates_local;
  size_t *index_duplicates_local;
  int     size_to_alloc_global=1024*1024;
  int     size_to_alloc_local =1024*1024;
  buffer_id             =(size_t *)SID_malloc(sizeof(size_t)*n_buffer_max);
  buffer_index_local    =(size_t *)SID_malloc(sizeof(size_t)*n_buffer_max);
  buffer_flag_local     =(char   *)SID_malloc(sizeof(char)  *n_buffer_max);
  id_duplicates_local   =(size_t *)SID_malloc(sizeof(size_t)*size_to_alloc_local);
  index_duplicates_local=(size_t *)SID_malloc(sizeof(size_t)*size_to_alloc_local);

  //... loop over all the ids on all ranks in buffer-sized batches
  SID_log("Checking for duplicates...",SID_LOG_OPEN|SID_LOG_TIMER);
  size_t i_particle;
  size_t j_particle;
  int    count;
  id_last=28991029248+1; // sufficient for Tiamat
  for(i_particle=0,count=0;i_particle<n_particles_all_in_groups;i_particle+=n_buffer){

    //Decide this buffer iteration's size
    n_buffer=MIN(n_buffer_max,n_particles_all_in_groups-i_particle);

    //Set the buffer to a default value smaller than (or equal to) the smallest possible data size
    for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
      buffer_id[i_buffer]         =0; // Min value is 0
      buffer_index_local[i_buffer]=n_particles_local+1;
      buffer_flag_local[i_buffer] =FALSE;
    }

    // Determine if any of the local data is needed for this buffer ...
    for(j_particle=0;j_particle<n_particles_local;j_particle++){
      index_rank_test_local=id_list_local_rank[j_particle]-i_particle;
      //... if so, set the appropriate buffer value
      if(index_rank_test_local>=0 && index_rank_test_local<n_buffer){
        buffer_id[index_rank_test_local]         =id_list_local[j_particle];
        buffer_index_local[index_rank_test_local]=j_particle;
        buffer_flag_local[index_rank_test_local] =TRUE; // This array lets ranks know which particles in the scan are local
      }
    }

    // Doing a global max on the buffer yields the needed buffer on all ranks
    SID_Allreduce(SID_IN_PLACE,buffer_id,n_buffer,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);

    // Check for duplicates
    for(i_buffer=0;i_buffer<n_buffer;i_buffer++){
      size_t id_test;
      switch(i_buffer){
         case 0:
           id_test=id_last; // Check across buffer boundaries
           break;
         default:
           id_test=buffer_id[i_buffer-1];
      }
      if(buffer_id[i_buffer]==id_test){
        // Make sure the first instance gets included in the list too
        if(count==0){
           n_duplicates_unique_global++;
           n_duplicates_global++;
           if(flag_last){
              id_duplicates_local[n_duplicates_local]   =id_last;
              index_duplicates_local[n_duplicates_local]=index_last;
              if(index_duplicates_local[n_duplicates_local]>=n_particles_local)
                 SID_trap_error("Invalid local duplicate index (1) (%lld; n_particles_local=%lld) in construct_duplicate_list().",
                                ERROR_LOGIC,index_duplicates_local[n_duplicates_local],n_particles_local);
              n_duplicates_local++;
              if(n_duplicates_local>=size_to_alloc_local){
                size_to_alloc_local   =(int)(1.5*(float)size_to_alloc_local);
                id_duplicates_local   =(size_t *)SID_realloc(id_duplicates_local,   sizeof(size_t)*size_to_alloc_local);
                index_duplicates_local=(size_t *)SID_realloc(index_duplicates_local,sizeof(size_t)*size_to_alloc_local);
              }
           }
        }
        // Construct the local arrays
        n_duplicates_global++;
        if(buffer_flag_local[i_buffer]){
           id_duplicates_local[n_duplicates_local]   =buffer_id[i_buffer];
           index_duplicates_local[n_duplicates_local]=buffer_index_local[i_buffer];
           if(index_duplicates_local[n_duplicates_local]>=n_particles_local)
              SID_trap_error("Invalid local duplicate index (2) (%lld; n_particles_local=%lld) in construct_duplicate_list().",
                             ERROR_LOGIC,index_duplicates_local[n_duplicates_local],n_particles_local);
           n_duplicates_local++;
           if(n_duplicates_local>=size_to_alloc_local){
             size_to_alloc_local   =(int)(1.5*(float)size_to_alloc_local);
             id_duplicates_local   =(size_t *)SID_realloc(id_duplicates_local,   sizeof(size_t)*size_to_alloc_local);
             index_duplicates_local=(size_t *)SID_realloc(index_duplicates_local,sizeof(size_t)*size_to_alloc_local);
           }
        }
        count++;
      }
      else
         count=0;
      id_last   =buffer_id[i_buffer];
      index_last=buffer_index_local[i_buffer];
      flag_last =buffer_flag_local[i_buffer];
    } // i_buffer
  }
  SID_log("No. of duplicates found = %d",SID_LOG_COMMENT,n_duplicates_global);
  SID_log("No. that are unique     = %d",SID_LOG_COMMENT,n_duplicates_unique_global);

  // Determine the group/subgroup indices of the local duplicates (if there is any)
  int    *index_duplicates_group_local   =NULL;
  int    *index_duplicates_subgroup_local=NULL;
  size_t *index_duplicates_local_index   =NULL;
  if(n_duplicates_local>0){
     
     // Fetch the halo catalog group/subgroup information. 
     int     n_groups_all;
     int     n_groups;
     int     n_subgroups_all;
     int     n_subgroups;
     int    *n_subgroups_group;
     size_t *offset_groups;
     size_t *offset_subgroup;
     int    *size_groups;
     int    *size_subgroups;
     n_groups_all     =((int    *)ADaPS_fetch(plist->data,"n_groups_all_%03d",            snapshot_number))[0];
     n_groups         =((int    *)ADaPS_fetch(plist->data,"n_groups_%03d",                snapshot_number))[0];
     n_subgroups_all  =((int    *)ADaPS_fetch(plist->data,"n_subgroups_all_%03d",         snapshot_number))[0];
     n_subgroups      =((int    *)ADaPS_fetch(plist->data,"n_subgroups_%03d",             snapshot_number))[0];
     n_subgroups_group= (int    *)ADaPS_fetch(plist->data,"n_subgroups_group_%03d",       snapshot_number);
     offset_groups    = (size_t *)ADaPS_fetch(plist->data,"particle_offset_group_%03d",   snapshot_number);
     offset_subgroup  = (size_t *)ADaPS_fetch(plist->data,"particle_offset_subgroup_%03d",snapshot_number);
     size_groups      = (int    *)ADaPS_fetch(plist->data,"n_particles_group_%03d",       snapshot_number);
     size_subgroups   = (int    *)ADaPS_fetch(plist->data,"n_particles_subgroup_%03d",    snapshot_number);
   
     // Sort the local duplicate index list
     sort(index_duplicates_local,(size_t)n_duplicates_local,&index_duplicates_local_index,SID_SIZE_T,SORT_LOCAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE); 
   
     // Determine which group/subgroup each duplicate belongs to
     int    i_group;
     int    i_subgroup;
     int    i_duplicate;
     size_t index_i;
     size_t index_j;
     index_duplicates_group_local   =(int *)SID_malloc(sizeof(int)*n_particles_local);
     index_duplicates_subgroup_local=(int *)SID_malloc(sizeof(int)*n_particles_local);
     for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++){
        index_duplicates_group_local[i_duplicate]   =-1;
        index_duplicates_subgroup_local[i_duplicate]=-1;
     }

     for(i_duplicate=0,i_group=0,i_subgroup=0;i_duplicate<n_duplicates_local;i_duplicate++){
        index_i=index_duplicates_local_index[i_duplicate];
        index_j=index_duplicates_local[index_i];
        // Identify group
        while(index_j>=(offset_groups[i_group]+(size_t)size_groups[i_group]) && i_group<(n_groups-1)) 
           i_group++;
        if(index_j>=offset_groups[i_group] && index_j<(offset_groups[i_group]+(size_t)size_groups[i_group]))
           index_duplicates_group_local[index_i]=i_group;
        // Identify subgroup
        while(index_j>=(offset_subgroup[i_subgroup]+(size_t)size_subgroups[i_subgroup]) && i_subgroup<(n_subgroups-1)) 
           i_subgroup++;
        if(index_j>=offset_subgroup[i_subgroup] && index_j<(offset_subgroup[i_subgroup]+(size_t)size_subgroups[i_subgroup]))
           index_duplicates_subgroup_local[index_i]=i_subgroup;
     }

     // Sanity check
     int n_bad=0;
     for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++){
        if(index_duplicates_group_local[i_duplicate]<0 && index_duplicates_subgroup_local[i_duplicate]<0)
           n_bad++;
     }
     SID_Allreduce(SID_IN_PLACE,&n_bad,1,SID_INT,SID_SUM,SID.COMM_WORLD);
     if(n_bad>0)
        SID_trap_error("%d duplicate(s) have not been sucessfully assigned to a group or subgroup.",ERROR_LOGIC,n_bad);
  }

  // Clean-up
  SID_free(SID_FARG buffer_id);
  SID_free(SID_FARG buffer_index_local);
  SID_free(SID_FARG buffer_flag_local);
  SID_free(SID_FARG id_list_local_rank);
  (*n_duplicates_global_return)            =n_duplicates_global;
  (*n_duplicates_local_return)             =n_duplicates_local;
  (*id_duplicates_local_return)            =id_duplicates_local;
  (*index_duplicates_local_return)         =index_duplicates_local;          // index in particle list
  (*index_duplicates_local_index_return)   =index_duplicates_local_index;    // sort indices of index in particle list
  (*index_duplicates_group_local_return)   =index_duplicates_group_local;    // index in group list
  (*index_duplicates_subgroup_local_return)=index_duplicates_subgroup_local; // index in subgroup list
  SID_log("Done.",SID_LOG_CLOSE);
}

void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              plist_info *plist,
                              size_t     *id_list,
                              size_t      n_particles_rank);
void read_gadget_binary_local(char       *filename_root_in,
                              int         snapshot_number,
                              plist_info *plist,
                              size_t     *id_list,
                              size_t      n_particles_rank){
  char    **pname;
  char     *name_initpositions;
  char      filename[MAX_FILENAME_LENGTH];
  char     *read_catalog;
  size_t    i,j,k,l,jj;
  size_t    n_particles_file;
  size_t    n_particles_kept;
  size_t    n_of_type_rank[N_GADGET_TYPE];
  size_t    n_of_type[N_GADGET_TYPE];
  int       n_type_used;
  int       n_particles_in_groups;
  size_t    n_particles_all;
  size_t    i_p,i_p_temp;
  size_t    id_min,id_max;
  double    mass_array[N_GADGET_TYPE];
  int       unused[GADGET_HEADER_SIZE];
  size_t    n_all[N_GADGET_TYPE];
  int       n_all_tmp[N_GADGET_TYPE];
  int       flag_metals;
  int       flag_ages;
  int       flag_entropyICs;
  double    expansion_factor;
  double    box_size;
  double    d_value;
  double    d1_value;
  double    d2_value;
  double    d3_value;
  float     f_value;
  int       i_value;
  int       i1_value;
  int       i2_value;
  int       i3_value;
  int       n_warning;
  GBPREAL  *x_array;
  GBPREAL  *y_array;
  GBPREAL  *z_array;
  size_t   *id_array;
  size_t   *id_list_offset;
  size_t   *id_list_in;
  size_t   *id_list_index;
  size_t    n_id_list;
  unsigned int       record_length_open;
  unsigned int       record_length_close;
  int       n_return;
  size_t    s_load;
  char     *keep;
  int       flag_keep_IDs      =TRUE;
  int       flag_multifile     =FALSE;
  int       flag_multimass     =FALSE;
  int       flag_file_type     =0;
  int       flag_filefound     =FALSE;
  int       flag_initpositions =FALSE;
  int       flag_LONGIDS       =FALSE;
  int       flag_read_marked   =FALSE;
  int       flag_read_catalog  =FALSE;
  int       flag_all_read_all  =FALSE;
  int       i_file;
  int       i_rank;
  int       n_files;
  int       n_files_in;
  double    x_scale;
  double    y_scale;
  double    z_scale;
  double    x_offset;
  double    y_offset;
  double    z_offset;
  double    x_period;
  double    y_period;
  double    z_period;
  int       flag_xscale;
  int       flag_yscale;
  int       flag_zscale;
  int       flag_xoffset;
  int       flag_yoffset;
  int       flag_zoffset;
  int       flag_xperiod;
  int       flag_yperiod;
  int       flag_zperiod;
  FILE     *fp;
  void     *buffer;
  size_t   *buffer_index=NULL;
  void     *buffer_i;
  size_t    id_search;
  size_t    id_value;
  size_t    n_buffer;
  int       n_non_zero;
  int       seed=1073743;
  int       read_mode;
  int       mark_mode;
  double    x_min_read;
  double    x_max_read;
  double    y_min_read;
  double    y_max_read;
  double    z_min_read;
  double    z_max_read;
  double    x_min_bcast;
  double    x_max_bcast;
  double    y_min_bcast;
  double    y_max_bcast;
  double    z_min_bcast;
  double    z_max_bcast;
  gadget_header_info header;

  flag_filefound=init_gadget_read(filename_root_in,snapshot_number,&flag_multifile,&flag_file_type,&header);

  // A file was found ... 
  SID_log("Reading GADGET binary file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in);
  if(flag_filefound){

    pname=plist->species;

    // Number of particles for each species in this file
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type[i]=(size_t)header.n_file[i];

    // Expansion factor (or time) 
    expansion_factor=header.time;
    ADaPS_store(&(plist->data),(void *)(&expansion_factor),"expansion_factor",ADaPS_SCALAR_DOUBLE);
    ADaPS_store(&(plist->data),(void *)(&(header.time)),   "time",            ADaPS_SCALAR_DOUBLE);

    // Redshift
    d_value=(double)header.redshift;
    ADaPS_store(&(plist->data),(void *)(&d_value),"redshift",ADaPS_SCALAR_DOUBLE);

    // Number of particles for each species in all files
    for(i=0;i<N_GADGET_TYPE;i++)
      n_all[i]=(size_t)header.n_all[i];

    // Number of files in this snapshot 
    ADaPS_store(&(plist->data),(void *)(&(header.n_files)),"n_files",ADaPS_SCALAR_INT);
    if(flag_multifile)
      n_files=header.n_files;

    // Cosmology

    // Omega_o
    d_value=(double)header.Omega_M;
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_M",ADaPS_SCALAR_DOUBLE);

    // Omega_Lambda 
    d_value=(double)header.Omega_Lambda;
    ADaPS_store(&(plist->data),(void *)(&d_value),"Omega_Lambda",ADaPS_SCALAR_DOUBLE);

    // Box size
    ADaPS_store(&(plist->data),(void *)(&box_size),"box_size",ADaPS_SCALAR_DOUBLE);

    // Species total particle counts and mass array
    for(i=0;i<N_GADGET_TYPE;i++){
      mass_array[i]=header.mass_array[i];
      if(header.n_all_high_word[i]>0)
        n_all[i]+=(((size_t)(header.n_all_high_word[i])) << 32);
    }

    // Count the total number of particles
    n_particles_all=0;
    for(i=0,n_non_zero=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        n_particles_all+=n_all[i];
        n_non_zero++;
      }
    }

    // Report the number of particles
    SID_log("%lld",SID_LOG_CONTINUE,n_particles_all);
    if(n_non_zero>0)
      SID_log(" (",SID_LOG_CONTINUE,n_particles_all);
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        if(i==n_non_zero-1){
          if(n_non_zero>1)
            SID_log("and %lld %s",SID_LOG_CONTINUE,n_all[i],pname[i]);
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_all[i],pname[i]);
        }
        else{
          if(n_non_zero>1)
            SID_log("%lld %s, ",SID_LOG_CONTINUE,n_all[i],pname[i]);
          else
            SID_log("%lld %s",SID_LOG_CONTINUE,n_all[i],pname[i]);
        }
      }
    }
    if(n_non_zero>0)
      SID_log(") particles...",SID_LOG_CONTINUE);
    else
      SID_log(" particles...",SID_LOG_CONTINUE);

    // Store mass array
    for(i=0;i<N_GADGET_TYPE;i++){
      if(n_all[i]>0){
        if(mass_array[i]!=0.){
          d_value=mass_array[i];
          ADaPS_store(&(plist->data),(void *)(&d_value),"mass_array_%s",ADaPS_SCALAR_DOUBLE,plist->species[i]);
        }
      }
      else
        mass_array[i]=0.;
    }

    // Fetch (and sort by id) local group particles
    SID_log("Sorting group particle ids...",SID_LOG_OPEN);
    merge_sort((void *)id_list,(size_t)n_particles_rank,&id_list_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
    for(i=0;i<N_GADGET_TYPE;i++)
      n_of_type_rank[i]=0;
    n_of_type_rank[GADGET_TYPE_DARK]=n_particles_rank;
    SID_log("Done.",SID_LOG_CLOSE);

    // Allocate data arrays
    id_array=(size_t  *)SID_malloc(sizeof(size_t) *(size_t)n_particles_rank);
    x_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*(size_t)n_particles_rank);
    y_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*(size_t)n_particles_rank);
    z_array =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*(size_t)n_particles_rank);

    // Read data
    for(i_file=0,n_particles_kept=0;i_file<n_files;i_file++){
      if(n_files>1)
        SID_log("Reading file %d of %d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file+1,n_files);
      else
        SID_log("Performing read...",SID_LOG_OPEN|SID_LOG_TIMER);

      // Set filename
      set_gadget_filename(filename_root_in,snapshot_number,i_file,flag_multifile,flag_file_type,filename);

      // Read header and set positions pointer
      FILE *fp_positions;
      if(SID.I_am_Master){
        fp_positions=fopen(filename,"r");
        fread(&record_length_open, sizeof(int),1,               fp_positions);
        fread(&header,             sizeof(gadget_header_info),1,fp_positions);
        fread(&record_length_close,sizeof(int),1,               fp_positions);
        if(record_length_open!=record_length_close)
          SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
        fread(&record_length_open,sizeof(int),1,fp_positions);
      }
      SID_Bcast(&header,sizeof(gadget_header_info),MASTER_RANK,SID.COMM_WORLD);
      for(i=0,n_particles_file=0;i<N_GADGET_TYPE;i++)
        n_particles_file+=(size_t)header.n_file[i];

      // Set IDs pointer
      FILE *fp_IDs;
      if(SID.I_am_Master){
         fp_IDs=fopen(filename,"r");
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         fseeko(fp_IDs,(off_t)record_length_open,SEEK_CUR);
         fread(&record_length_close,sizeof(int),1,fp_IDs);
         if(record_length_open!=record_length_close)
           SID_log_warning("Problem with GADGET record size (close of header)",ERROR_LOGIC);
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         fseeko(fp_IDs,(off_t)record_length_open,SEEK_CUR);
         fread(&record_length_close,sizeof(int),1,fp_IDs);
         if(record_length_open!=record_length_close)
           SID_log_warning("Problem with GADGET record size (close of positions)",ERROR_LOGIC);
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         fseeko(fp_IDs,(off_t)record_length_open,SEEK_CUR);
         fread(&record_length_close,sizeof(int),1,fp_IDs);
         if(record_length_open!=record_length_close)
           SID_log_warning("Problem with GADGET record size (close of velocities)",ERROR_LOGIC);
         fread(&record_length_open,sizeof(int),1,fp_IDs);
         if((size_t)record_length_open/n_particles_file==sizeof(long long)){
            SID_log("(long long IDs)...",SID_LOG_CONTINUE);
            flag_LONGIDS=TRUE;
         }
         else if((size_t)record_length_open/n_particles_file==sizeof(int)){
            SID_log("(int IDs)...",SID_LOG_CONTINUE);
            flag_LONGIDS=FALSE;
         }
         else
            SID_trap_error("IDs record length (%d) does not set a sensible id byte size for the number of particles given in the header (%s)",ERROR_LOGIC,record_length_open,n_particles_file);
      }
      SID_Bcast(&flag_LONGIDS,sizeof(int),MASTER_RANK,SID.COMM_WORLD);

      // Allocate buffers
      int      n_buffer_max=MIN(n_particles_file,4*1024*1024);
      GBPREAL *buffer_positions ;
      size_t  *buffer_IDs       ;
      buffer_positions =(GBPREAL *)SID_malloc(3*n_buffer_max*sizeof(GBPREAL));
      buffer_IDs       =(size_t  *)SID_malloc(  n_buffer_max*sizeof(size_t));

      // Perform read
      int     i_buffer;
      int     n_buffer_left;
      int     n_buffer;
      int     i_species;
      int     i_particle;
      int     j_particle;
      int     i_list;
      size_t *buffer_index=NULL;
      for(i_species=0,i_buffer=n_buffer_max,i_list=0,n_buffer_left=n_particles_file;i_species<N_GADGET_TYPE;i_species++){
         if(header.n_file[i_species]>0){
            for(i_particle=0;i_particle<header.n_file[i_species];i_particle++,i_buffer++){
               // Perform a buffered read
               if(i_buffer>=n_buffer_max){
                  n_buffer=MIN(n_buffer_max,n_buffer_left);
                  if(SID.I_am_Master){
                     fread(buffer_positions, sizeof(GBPREAL),3*n_buffer,fp_positions);
                     if(!flag_LONGIDS){
                        int *buffer_IDs_int=(int *)buffer_IDs;
                        fread(buffer_IDs_int,sizeof(int),n_buffer,fp_IDs);
                        for(j_particle=n_buffer-1;j_particle>=0;j_particle--)
                           buffer_IDs[j_particle]=(size_t)buffer_IDs_int[j_particle];
                     }
                     else
                        fread(buffer_IDs,sizeof(size_t),n_buffer,fp_IDs);
                  }
                  SID_Bcast(buffer_positions, 3*n_buffer*sizeof(GBPREAL),MASTER_RANK,SID.COMM_WORLD);
                  SID_Bcast(buffer_IDs,         n_buffer*sizeof(size_t), MASTER_RANK,SID.COMM_WORLD);
                  SID_free(SID_FARG buffer_index);
                  merge_sort(buffer_IDs,(size_t)n_buffer,&buffer_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
                  n_buffer_left-=n_buffer;
                  i_buffer=0;
                  i_list  =0;
               }
               if(i_list<n_particles_rank){
                  // Move to the next particle
                  size_t id_test;
                  size_t j_buffer=buffer_index[i_buffer];
                  id_test=buffer_IDs[j_buffer];
                  while(id_list[id_list_index[i_list]]<id_test && i_list<(n_particles_rank-1))
                    i_list++;

                  // If we've found a particle, add it to the arrays
                  if(i_list<n_particles_rank){
                    while(id_list[id_list_index[i_list]]==id_test){ // The while takes care of duplicates
                      GBPREAL f1_value;
                      GBPREAL f2_value;
                      GBPREAL f3_value;
                      f1_value=((GBPREAL *)buffer_positions)[3*j_buffer+0];
                      f2_value=((GBPREAL *)buffer_positions)[3*j_buffer+1];
                      f3_value=((GBPREAL *)buffer_positions)[3*j_buffer+2];
                      id_array[id_list_index[i_list]]=id_test;
                      x_array[id_list_index[i_list]] =f1_value;
                      y_array[id_list_index[i_list]] =f2_value;
                      z_array[id_list_index[i_list]] =f3_value;
                      n_particles_kept++;
                      i_list++;
                      if(i_list>=n_particles_rank) break;
                    }
                  }
               }
            } // i_particle
         } // if n>0
      } // i_species

      // Clean-up
      SID_free(SID_FARG buffer_positions);
      SID_free(SID_FARG buffer_IDs);
      SID_free(SID_FARG buffer_index);
      if(SID.I_am_Master){
         fclose(fp_positions);
         fclose(fp_IDs);
      }
      SID_log("Done.",SID_LOG_CLOSE);
    } // i_file
    SID_free(SID_FARG id_list_index);

    // Check that the right number of particles have been read
    if(n_particles_kept!=n_particles_rank)
      SID_log_warning("Rank %d did not receive the right number of particles (ie. %d!=%d)",ERROR_LOGIC,SID.My_rank,n_particles_kept,n_particles_rank);
    SID_Allreduce(&n_particles_rank,&n_particles_all,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);

    // Store everything in the data structure...
    //   ... particle counts ...
    ADaPS_store(&(plist->data),(void *)(&n_particles_rank),"n_list",    ADaPS_SCALAR_SIZE_T);
    ADaPS_store(&(plist->data),(void *)(&n_particles_all), "n_all_list",ADaPS_SCALAR_SIZE_T);

    //   ... positions ...
    ADaPS_store(&(plist->data),(void *)x_array,"x_list",ADaPS_DEFAULT);
    ADaPS_store(&(plist->data),(void *)y_array,"y_list",ADaPS_DEFAULT);
    ADaPS_store(&(plist->data),(void *)z_array,"z_list",ADaPS_DEFAULT);

    //  ... ids ...
    ADaPS_store(&(plist->data),(void *)id_array,"id_list",ADaPS_DEFAULT);

    SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not find file with root {%s}",ERROR_IO_OPEN,filename_root_in);
}

int main(int argc, char *argv[]){
  plist_info  plist;
  char        filename_halos_in_root[256];
  char        filename_halos_out_root[256];
  char        filename_snapshot_root[256];
  char        filename_snapshot[256];
  char        filename_output_properties_dir[256];
  char        filename_output_properties[256];
  char        filename_output_properties_temp[256];
  char        group_text_prefix[4];
  int         n_groups_process;
  int         n_groups_local;
  int         n_groups_all;
  int         i_rank;
  int         i_group;
  int         i_file_lo_in;
  int         i_file_lo;
  int         i_file_hi_in;
  int         i_file_hi;
  int         i_file;
  int         i_file_skip;
  int         i_particle;
  int         j_particle;
  int         i_process;
  int         n_particles;
  int         n_particles_max;
  size_t      n_particles_cumulative;
  size_t     *ids_particles;
  size_t     *ids_particles_index;
  size_t     *ids_groups;
  int        *n_particles_groups_process;
  int        *n_particles_groups;
  int        *n_particles_subgroups;
  int        *group_offset;
  size_t      n_particles_in_groups;
  size_t     *ids_snapshot;
  size_t     *ids_sort_index;
  size_t     *ids_snapshot_sort_index;
  size_t     *ids_groups_sort_index;
  size_t      n_particles_snapshot;
  int         i_value;
  double      Omega_M;
  double      Omega_b;
  double      Omega_Lambda;
  double      f_gas;
  double      Omega_k;
  double      sigma_8;
  double      n_spec;
  double      redshift;
  double      box_size;
  double      particle_mass;
  
  int         r_val;
  struct      stat file_stats;
  size_t      n_bytes;
  size_t      n_bytes_buffer;
  cosmo_info *cosmo;
  int                   n_temp;
  int                   n_truncated;
  int                   largest_truncated;
  int                   largest_truncated_local;
  char                 *filename_number;

  SID_init(&argc,&argv,NULL);

  // Fetch user inputs
  strcpy(filename_halos_in_root, argv[1]);
  strcpy(filename_snapshot_root, argv[2]);
  strcpy(filename_halos_out_root,argv[3]);
  box_size    =atof(argv[4]);
  i_file_lo_in=atoi(argv[5]);
  i_file_hi_in=atoi(argv[6]);
  i_file_skip =atoi(argv[7]);

  if(i_file_lo_in<i_file_hi_in){
     i_file_lo=i_file_lo_in;
     i_file_hi=i_file_hi_in;
  }
  else{
     i_file_lo=i_file_hi_in;
     i_file_hi=i_file_lo_in;
  }

  SID_log("Removing duplicates from halo files for snapshots #%d->#%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file_lo,i_file_hi);
  for(i_file=i_file_hi;i_file>=i_file_lo;i_file-=i_file_skip){
    SID_log("Processing file #%03d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,1);

    // Read group info from the halo catalogs
    int    *n_subgroups_group;
    size_t *offset_groups;
    size_t *offset_subgroups;
    int    *size_groups;
    int    *size_subgroups;
    int     n_subgroups_all;
    int     n_subgroups_local;
    init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
    filename_number=(char *)SID_malloc(sizeof(char)*10);
    sprintf(filename_number,"%03d", i_file);
    ADaPS_store(&(plist.data),(void *)filename_number,"read_catalog",ADaPS_DEFAULT);
    read_groups(filename_halos_in_root,i_file,READ_GROUPS_ALL,&plist,filename_number);

    // Fetch results of read_groups()
    n_groups_all     =((int    *)ADaPS_fetch(plist.data,"n_groups_all_%03d",             i_file))[0];
    n_groups_local   =((int    *)ADaPS_fetch(plist.data,"n_groups_%03d",                 i_file))[0];
    n_subgroups_all  =((int    *)ADaPS_fetch(plist.data,"n_subgroups_all_%03d",          i_file))[0];
    n_subgroups_local=((int    *)ADaPS_fetch(plist.data,"n_subgroups_%03d",              i_file))[0];
    n_subgroups_group= (int    *)ADaPS_fetch(plist.data,"n_subgroups_group_%03d",        i_file);
    n_subgroups_group= (int    *)ADaPS_fetch(plist.data,"n_subgroups_group_%03d",        i_file);
    offset_groups    = (size_t *)ADaPS_fetch(plist.data,"particle_offset_group_%03d",    i_file);
    offset_subgroups = (size_t *)ADaPS_fetch(plist.data,"particle_offset_subgroup_%03d", i_file);
    size_groups      = (int    *)ADaPS_fetch(plist.data,"n_particles_group_%03d",        i_file);
    size_subgroups   = (int    *)ADaPS_fetch(plist.data,"n_particles_subgroup_%03d",     i_file);

    // Determine the offsets necessary to place local offsets/indices to global values
    SID_log("Determine global offsets of the halo catalog...",SID_LOG_OPEN|SID_LOG_TIMER);
    int group_index_local_offset        =0;
    int subgroup_index_local_offset     =0;
    int group_index_local_offset_comm   =0;
    int subgroup_index_local_offset_comm=0;
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          group_index_local_offset        =group_index_local_offset_comm;
          subgroup_index_local_offset     =subgroup_index_local_offset_comm;
          group_index_local_offset_comm   =group_index_local_offset   +n_groups_local;
          subgroup_index_local_offset_comm=subgroup_index_local_offset+n_groups_local;
       }
       SID_Bcast(&group_index_local_offset_comm,   sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(&subgroup_index_local_offset_comm,sizeof(int),i_rank,SID.COMM_WORLD);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // For testing purposes: read a list of duplicates from the ./remove_duplicates_testing.dat
    //    file and inject them into the catalog.  Does nothing if the file is not present.
    //    The file should have 4 columns of integers -- one line per artificial duplicate --
    //    with the first two columns specifying the group number and group particle index 
    //    for the particle to be duplicated, followed by the group number and group particle
    //    index of where it will be duplicated to.
    int     flag_add_duplicates;
    int     n_duplicates_added;
    int     n_duplicates_added_local;
    size_t *id_duplicate_added;
    int    *group_duplicate_added_to;
    int    *subgroup_duplicate_added_to;
    int    *index_duplicate_added_to;
    read_duplicates_file_local("./remove_duplicates_testing.dat",filename_number,
                               &plist,
                               size_groups,
                               offset_groups,
                               size_subgroups,
                               offset_subgroups,
                               n_subgroups_group,
                               n_groups_local,
                               n_subgroups_local,
                               group_index_local_offset,
                               subgroup_index_local_offset,
                               &flag_add_duplicates,
                               &n_duplicates_added,
                               &n_duplicates_added_local,
                               &id_duplicate_added,
                               &group_duplicate_added_to,
                               &subgroup_duplicate_added_to,
                               &index_duplicate_added_to);
    if(flag_add_duplicates){
       inject_duplicates_local(&plist,filename_number,
                               size_groups,
                               offset_groups,
                               size_subgroups,
                               offset_subgroups,
                               n_groups_local,
                               n_subgroups_local,
                               group_index_local_offset,
                               subgroup_index_local_offset,
                               n_duplicates_added,
                               n_duplicates_added_local,
                               id_duplicate_added,
                               group_duplicate_added_to,
                               subgroup_duplicate_added_to,
                               index_duplicate_added_to);
       SID_free(SID_FARG id_duplicate_added);
       SID_free(SID_FARG group_duplicate_added_to);
       SID_free(SID_FARG subgroup_duplicate_added_to);
       SID_free(SID_FARG index_duplicate_added_to);
    }

    // These things could have been changed during duplicate injection, so get them now
    size_t *id_list_local;
    size_t  n_particles_all;
    size_t  n_particles_local;
    n_particles_all  =((size_t *)ADaPS_fetch(plist.data,"n_particles_all_%03d",i_file))[0];
    n_particles_local=((size_t *)ADaPS_fetch(plist.data,"n_particles_%03d",    i_file))[0];
    id_list_local    = (size_t *)ADaPS_fetch(plist.data,"particle_ids_%03d",   i_file);

    // If there's any groups to analyze ...
    size_t *id_duplicates_local            =NULL;
    size_t *index_duplicates_local         =NULL;
    size_t *index_duplicates_local_index   =NULL;
    int    *index_duplicates_group_local   =NULL;
    int    *index_duplicates_subgroup_local=NULL;
    int     n_duplicates_global            =0;
    int     n_duplicates_local             =0;
    if(n_groups_all>0){

       // Determine duplicates
       construct_duplicate_list_local(&plist,
                                      i_file,
                                      id_list_local,
                                      n_particles_local,
                                      n_particles_all,
                                      &id_duplicates_local,             // sorted
                                      &index_duplicates_local,          // particle index of duplicate
                                      &index_duplicates_local_index,
                                      &index_duplicates_group_local,    // group index of duplicate
                                      &index_duplicates_subgroup_local, // subgroup index of duplicate
                                      &n_duplicates_global,
                                      &n_duplicates_local);

       // Create a master list of all IDs we will grab positions for from the snapshot ...
       SID_log("Creating master list of particles including MBPs...",SID_LOG_OPEN);

       // ... first: copy the duplicates to the start of the list ...
       size_t *id_list_master;
       int     n_id_list_master;
       int    *locator_MBP;
       n_id_list_master=2*n_duplicates_local; // One for the duplicate, one for the MBP it will be tested against
       id_list_master  =(size_t *)SID_malloc(sizeof(size_t)*n_id_list_master);
       memcpy(id_list_master,id_duplicates_local,sizeof(size_t)*n_duplicates_local);

       // ... second: add the local MBPs (subgroups first, then groups)
       int i_list;
       int i_group;
       int i_subgroup;
       int i_duplicate;
       for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++){
          if(index_duplicates_subgroup_local[i_duplicate]>=0 && index_duplicates_subgroup_local[i_duplicate]<n_subgroups_local)
             id_list_master[n_duplicates_local+i_duplicate]=id_list_local[offset_subgroups[index_duplicates_subgroup_local[i_duplicate]]];
          else if(index_duplicates_group_local[i_duplicate]>=0 && index_duplicates_group_local[i_duplicate]<n_groups_local)
             id_list_master[n_duplicates_local+i_duplicate]=id_list_local[offset_groups[index_duplicates_group_local[i_duplicate]]];
          else
             SID_trap_error("Neither a valid group nor valid subgroup has been specified for a duplicate particle.",ERROR_LOGIC);
       }
       SID_log("Done.",SID_LOG_CLOSE);
       
       // Read the needed particle positions from the snapshot
       read_gadget_binary_local(filename_snapshot_root,i_file,&plist,id_list_master,(size_t)n_id_list_master);

       // Compute each duplicate's separation from its MBP
       SID_log("Computing displacements from MBPs...",SID_LOG_OPEN|SID_LOG_TIMER);
       GBPREAL *x_array;
       GBPREAL *y_array;
       GBPREAL *z_array;
       GBPREAL *r2_duplicates_local;
       GBPREAL  dx,dy,dz;
       r2_duplicates_local=(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_duplicates_local);
       x_array            =(GBPREAL *)ADaPS_fetch(plist.data,"x_list");
       y_array            =(GBPREAL *)ADaPS_fetch(plist.data,"y_list");
       z_array            =(GBPREAL *)ADaPS_fetch(plist.data,"z_list");
       for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++){
          dx=d_periodic(x_array[i_duplicate]-x_array[n_duplicates_local+i_duplicate],box_size);
          dy=d_periodic(y_array[i_duplicate]-y_array[n_duplicates_local+i_duplicate],box_size);
          dz=d_periodic(z_array[i_duplicate]-z_array[n_duplicates_local+i_duplicate],box_size);
          r2_duplicates_local[i_duplicate]=dx*dx+dy*dy+dz*dz;
       }
       SID_log("Done.",SID_LOG_CLOSE);

       // Sort the IDs of the local duplicates
       size_t *id_duplicates_local_index=NULL;
       sort(id_duplicates_local,n_duplicates_local,&id_duplicates_local_index,SID_SIZE_T,SORT_LOCAL,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

       // Convert group/subgroup offsets to global references (instead of their current local reference)
       SID_log("Computing global particle offsets...",SID_LOG_OPEN);
       size_t offset_to_global_local=0;
       size_t offset_to_global_bcast=0;
       for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(SID.My_rank==i_rank){
             offset_to_global_local=offset_to_global_bcast;
             offset_to_global_bcast=offset_to_global_local+n_particles_local;
          }
          SID_Bcast(&offset_to_global_bcast,sizeof(size_t),i_rank,SID.COMM_WORLD);
       }
       SID_log("Done.",SID_LOG_CLOSE);

       // Decide which duplicates to keep and which to remove
       SID_log("Choose which duplicates to keep...",SID_LOG_OPEN|SID_LOG_TIMER);
       int      n_duplicates_exchange;
       size_t  *id_duplicates_exchange;
       size_t  *id_duplicates_exchange_index;
       GBPREAL *r2_duplicates_exchange;
       size_t   offset_to_global_exchange;
       char    *flag_keep_local;
       int      n_exchange_max;
       SID_Allreduce(&n_duplicates_local,&n_exchange_max,1,SID_INT,SID_MAX,SID.COMM_WORLD);
       if(SID.n_proc>1){
          id_duplicates_exchange      =(size_t  *)SID_malloc(sizeof(size_t) *n_exchange_max);
          id_duplicates_exchange_index=(size_t  *)SID_malloc(sizeof(size_t) *n_exchange_max);
          r2_duplicates_exchange      =(GBPREAL *)SID_malloc(sizeof(GBPREAL)*n_exchange_max);
       }
       flag_keep_local=(char *)SID_malloc(sizeof(char)*n_duplicates_local);
       for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++)
          flag_keep_local[i_duplicate]=TRUE;
       for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          // Perform ring exchanges
          if(i_rank==0){
             n_duplicates_exchange    =n_duplicates_local;
             offset_to_global_exchange=offset_to_global_local;
             if(SID.n_proc>1){
                memcpy(id_duplicates_exchange,      id_duplicates_local,      sizeof(size_t) *n_duplicates_local);
                memcpy(id_duplicates_exchange_index,id_duplicates_local_index,sizeof(size_t) *n_duplicates_local);
                memcpy(r2_duplicates_exchange,      r2_duplicates_local,      sizeof(GBPREAL)*n_duplicates_local);
             }
             else{
                id_duplicates_exchange      =id_duplicates_local;
                id_duplicates_exchange_index=id_duplicates_local_index;
                r2_duplicates_exchange      =r2_duplicates_local;
             }
          }
          else{
             exchange_ring_buffer(&n_duplicates_local,
                                  sizeof(int),
                                  1,
                                  &n_duplicates_exchange,
                                  NULL,
                                  i_rank);
             exchange_ring_buffer(&offset_to_global_local,
                                  sizeof(size_t),
                                  1,
                                  &offset_to_global_exchange,
                                  NULL,
                                  i_rank);
             exchange_ring_buffer(id_duplicates_local,
                                  sizeof(size_t),
                                  (size_t)n_duplicates_local,
                                  id_duplicates_exchange,
                                  NULL,
                                  i_rank);
             exchange_ring_buffer(id_duplicates_local_index,
                                  sizeof(size_t),
                                  (size_t)n_duplicates_local,
                                  id_duplicates_exchange_index,
                                  NULL,
                                  i_rank);
             exchange_ring_buffer(r2_duplicates_local,
                                  sizeof(GBPREAL),
                                  (size_t)n_duplicates_local,
                                  r2_duplicates_exchange,
                                  NULL,
                                  i_rank);
          }

          // Flag particles to be discarded
          int j_duplicate;
          for(i_duplicate=0;i_duplicate<n_duplicates_local && n_duplicates_exchange>0;i_duplicate++){
             j_duplicate=find_index(id_duplicates_exchange,
                                    id_duplicates_local[id_duplicates_local_index[i_duplicate]],
                                    n_duplicates_exchange,
                                    id_duplicates_exchange_index);
             while(id_duplicates_exchange[id_duplicates_exchange_index[j_duplicate]]<id_duplicates_local[id_duplicates_local_index[i_duplicate]]){
                j_duplicate++;
                if(j_duplicate<(n_duplicates_exchange-1))
                   continue;
                else
                   break;
             }
             while(id_duplicates_exchange[id_duplicates_exchange_index[j_duplicate]]==id_duplicates_local[id_duplicates_local_index[i_duplicate]]){
                // Prevent checks-against-self
                if(i_rank==0 && i_duplicate==j_duplicate){
                   j_duplicate++;
                   if(j_duplicate<n_duplicates_exchange)
                      continue;
                   else
                      break;
                }
                // Simplest case of there being a better particle
                if(r2_duplicates_exchange[id_duplicates_exchange_index[j_duplicate]]<r2_duplicates_local[id_duplicates_local_index[i_duplicate]])
                   flag_keep_local[id_duplicates_local_index[i_duplicate]]=FALSE;
                // More complicated case of there being an equally good particle (somehow)
                else if(r2_duplicates_exchange[id_duplicates_exchange_index[j_duplicate]]==r2_duplicates_local[id_duplicates_local_index[i_duplicate]]){
                  // Keep the duplicate instance that comes first in the particle list if displacements are identical
                  int index_global_duplicate_local;
                  int index_global_duplicate_exchange;
                  index_global_duplicate_local   =id_duplicates_local_index[i_duplicate]   +offset_to_global_local;
                  index_global_duplicate_exchange=id_duplicates_exchange_index[j_duplicate]+offset_to_global_exchange;
                  // Even if checks-against-self get here, they will go no further so long as we use '<' here and not '<='
                  if(index_global_duplicate_local<index_global_duplicate_exchange)
                     flag_keep_local[id_duplicates_local_index[i_duplicate]]=FALSE;
                }
                j_duplicate++;
                if(j_duplicate<n_duplicates_exchange)
                   continue;
                else
                   break;
             }
          }
       }
       if(SID.n_proc>1){
          SID_free(SID_FARG id_duplicates_exchange);
          SID_free(SID_FARG id_duplicates_exchange_index);
          SID_free(SID_FARG r2_duplicates_exchange);
       }
       SID_log("Done.",SID_LOG_CLOSE);

       // Write duplicates file
       write_duplicates_local(filename_halos_out_root,
                              filename_number,
                              n_duplicates_local,
                              id_duplicates_local,
                              index_duplicates_group_local,
                              index_duplicates_subgroup_local,
                              r2_duplicates_local,
                              flag_keep_local,
                              group_index_local_offset,
                              subgroup_index_local_offset);

       // Loop over the duplicates, removing the unwanted ones from the halo catalogs
       SID_log("Remove duplicates from catalog...",SID_LOG_OPEN|SID_LOG_TIMER);
       int n_removed      =0;
       int n_removed_local=0;
       for(i_particle=0,j_particle=0,i_duplicate=0,n_removed_local=0;
           i_particle<n_particles_local && i_duplicate<n_duplicates_local;i_duplicate++){
          // Move to the next duplicate
          size_t index_i;
          size_t index_j;
          index_i=index_duplicates_local_index[i_duplicate];
          index_j=index_duplicates_local[index_i]; // The duplicate indices shift down as we remove particles
          while(i_particle<index_j && i_particle<(n_particles_local-1))
             id_list_local[j_particle++]=id_list_local[i_particle++];
          // Sanity check
          if(id_list_local[i_particle]!=id_duplicates_local[index_i])
             SID_trap_error("Could not find a duplicate in the ID list.",ERROR_LOGIC);
          // Keep this duplicate ...
          if(flag_keep_local[index_i])
             id_list_local[j_particle++]=id_list_local[i_particle++];
          // ... or remove it.
          else{
             i_particle++;
             n_removed_local++;
          }
       }
       // Make sure the list is changed right through to the end
       while(i_particle<n_particles_local)
          id_list_local[j_particle++]=id_list_local[i_particle++];

       for(i_duplicate=0;i_duplicate<n_duplicates_local;i_duplicate++){
          if(!(flag_keep_local[i_duplicate])){
             if(index_duplicates_group_local[i_duplicate]>=0){
                size_groups[index_duplicates_group_local[i_duplicate]]--;
                for(i_group=index_duplicates_group_local[i_duplicate]+1;i_group<n_groups_local;i_group++)
                   offset_groups[i_group]--;
             }
             if(index_duplicates_subgroup_local[i_duplicate]>=0){
                size_subgroups[index_duplicates_subgroup_local[i_duplicate]]--;
                for(i_subgroup=index_duplicates_subgroup_local[i_duplicate]+1;i_subgroup<n_subgroups_local;i_subgroup++)
                   offset_subgroups[i_subgroup]--;
             }
          }
       }

       SID_free(SID_FARG r2_duplicates_local);
       SID_Allreduce(&n_removed_local,&n_removed,1,SID_INT,SID_SUM,SID.COMM_WORLD);
       n_particles_local-=n_removed_local;
       n_particles_all  -=n_removed;
       SID_log("No. of duplicates removed = %d",SID_LOG_COMMENT,n_removed);
       SID_log("Done.",SID_LOG_CLOSE);
    }

    // Convert group/subgroup offsets to global references (instead of their current local reference)
    SID_log("Determine global offsets of the halo catalog...",SID_LOG_OPEN|SID_LOG_TIMER);
    size_t offset_to_global=0;
    int    i_subgroup;
    int    j_subgroup;
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          for(i_group=0;i_group<n_groups_local;i_group++)
             offset_groups[i_group]+=offset_to_global;
          for(i_subgroup=0;i_subgroup<n_subgroups_local;i_subgroup++)
             offset_subgroups[i_subgroup]+=offset_to_global;
          offset_to_global+=n_particles_local;
       }
       SID_Bcast(&offset_to_global,sizeof(size_t),i_rank,SID.COMM_WORLD);
    }
    SID_log("Done.",SID_LOG_CLOSE);

    // Remove groups and subgroups that no longer have any particles
    SID_log("Remove zero-sized groups/subgroups...",SID_LOG_OPEN|SID_LOG_TIMER);
    int j_group;
    int n_groups_removed;
    int n_subgroups_removed;
    int n_subgroups_removed_i;
    for(i_group=0,j_group=0,i_subgroup=0,j_subgroup=0,n_subgroups_removed=0,n_groups_removed=0;i_group<n_groups_local;i_group++){
       int k_subgroup;
       for(k_subgroup=0,n_subgroups_removed_i=0;k_subgroup<n_subgroups_group[i_group];i_subgroup++,k_subgroup++){
          if(size_subgroups[i_subgroup]>0){
             offset_subgroups[j_subgroup]=offset_subgroups[i_subgroup];
             size_subgroups[j_subgroup]  =size_subgroups[i_subgroup];
             j_subgroup++;
          }
          else
             n_subgroups_removed_i++;
       }
       if(size_groups[i_group]>0){
          offset_groups[j_group]     =offset_groups[i_group];
          size_groups[j_group]       =size_groups[i_group];
          n_subgroups_group[j_group]=n_subgroups_group[i_group]-n_subgroups_removed_i;
          j_group++;
       }
       else
          n_groups_removed++;
       n_subgroups_removed+=n_subgroups_removed_i;
    }

    // Recompute the global group/subgroup counts
    n_groups_local   =j_group;
    n_subgroups_local=j_subgroup;
    calc_sum_global(&n_groups_local,   &n_groups_all,   1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    calc_sum_global(&n_subgroups_local,&n_subgroups_all,1,SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
    SID_log("No. of groups    removed = %d",SID_LOG_COMMENT,n_groups_removed);
    SID_log("No. of subgroups removed = %d",SID_LOG_COMMENT,n_subgroups_removed);
    SID_log("Done.",SID_LOG_CLOSE);

    // Allocate a buffer for writing
    int     n_buffer_alloc;
    void   *buffer;
    int    *buffer_int;
    size_t *buffer_size_t;
    n_buffer_alloc=n_particles_local;
    SID_Allreduce(SID_IN_PLACE,&n_buffer_alloc,1,SID_INT,SID_MAX,SID.COMM_WORLD);
    buffer       =SID_malloc(sizeof(size_t)*n_buffer_alloc);
    buffer_int   =(int    *)buffer;
    buffer_size_t=(size_t *)buffer;

    // Read the ID and offset byte-sizes from the input catalog
    FILE *fp_in_groups;
    FILE *fp_in_subgroups;
    FILE *fp_in_particles;
    char  filename_in_groups[MAX_FILENAME_LENGTH];
    char  filename_in_subgroups[MAX_FILENAME_LENGTH];
    char  filename_in_particles[MAX_FILENAME_LENGTH];
    int   n_byte_offsets_groups;
    int   n_byte_offsets_subgroups;
    int   n_byte_ids;
    sprintf(filename_in_groups,   "%s_%s.catalog_groups",   filename_halos_in_root,filename_number);
    sprintf(filename_in_subgroups,"%s_%s.catalog_subgroups",filename_halos_in_root,filename_number);
    sprintf(filename_in_particles,"%s_%s.catalog_particles",filename_halos_in_root,filename_number);
    if(SID.I_am_Master){
       int n_groups_in;
       int n_subgroups_in;
       fp_in_groups   =fopen(filename_in_groups,   "r");
       fp_in_subgroups=fopen(filename_in_subgroups,"r");
       fp_in_particles=fopen(filename_in_particles,"r");
       fread(&n_groups_in,             sizeof(int),1,fp_in_groups);
       fread(&n_byte_offsets_groups,   sizeof(int),1,fp_in_groups);
       fread(&n_subgroups_in,          sizeof(int),1,fp_in_subgroups);
       fread(&n_byte_offsets_subgroups,sizeof(int),1,fp_in_subgroups);
       fread(&n_byte_ids,              sizeof(int),1,fp_in_particles);
       fclose(fp_in_groups);
       fclose(fp_in_subgroups);
       fclose(fp_in_particles);
    }
    SID_Bcast(&n_byte_offsets_groups,   sizeof(int),MASTER_RANK,SID.COMM_WORLD);
    SID_Bcast(&n_byte_offsets_subgroups,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
    SID_Bcast(&n_byte_ids,              sizeof(int),MASTER_RANK,SID.COMM_WORLD);
    SID_log("Halo catalog byte-sizes:", SID_LOG_OPEN);
    SID_log("particle IDs    = %d byte",SID_LOG_COMMENT,n_byte_ids);
    SID_log("group    offset = %d byte",SID_LOG_COMMENT,n_byte_offsets_groups);
    SID_log("subgroup offset = %d byte",SID_LOG_COMMENT,n_byte_offsets_subgroups);
    SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);

    // Write results
    SID_log("Writing new halo catalog files...",SID_LOG_OPEN|SID_LOG_TIMER);
    int   i_rank;
    int   buffer_size;
    FILE *fp_out_groups;
    FILE *fp_out_subgroups;
    FILE *fp_out_particles;
    char  filename_out_groups[MAX_FILENAME_LENGTH];
    char  filename_out_subgroups[MAX_FILENAME_LENGTH];
    char  filename_out_particles[MAX_FILENAME_LENGTH];
    sprintf(filename_out_groups,   "%s_%s.catalog_groups",   filename_halos_out_root,filename_number);
    sprintf(filename_out_subgroups,"%s_%s.catalog_subgroups",filename_halos_out_root,filename_number);
    sprintf(filename_out_particles,"%s_%s.catalog_particles",filename_halos_out_root,filename_number);

    // Open files and write header
    if(SID.I_am_Master){
       // Open files
       fp_out_groups   =fopen(filename_out_groups,   "w");
       fp_out_subgroups=fopen(filename_out_subgroups,"w");
       fp_out_particles=fopen(filename_out_particles,"w");
       // Write headers
       fwrite(&n_groups_all,            sizeof(int),   1,fp_out_groups);
       fwrite(&n_byte_offsets_groups,   sizeof(int),   1,fp_out_groups);
       fwrite(&n_subgroups_all,         sizeof(int),   1,fp_out_subgroups);
       fwrite(&n_byte_offsets_subgroups,sizeof(int),   1,fp_out_subgroups);
       fwrite(&n_byte_ids,              sizeof(int),   1,fp_out_particles);
       if(n_byte_ids==sizeof(unsigned int)){
          int n_particles_head=(int)n_particles_all;
          fwrite(&n_particles_head,sizeof(int),1,fp_out_particles);
       }
       else{
          int64_t n_particles_head=(int64_t)n_particles_all;
          fwrite(&n_particles_head,sizeof(int64_t),1,fp_out_particles);
       }
    }

    // Write group sizes
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          buffer_size=(int)(sizeof(int)*n_groups_local);
          memcpy(buffer,size_groups,buffer_size);
       }
       SID_Bcast(&buffer_size,sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,      buffer_size,i_rank,SID.COMM_WORLD); 
       if(SID.I_am_Master){
          int n_write;
          n_write=fwrite(buffer,1,buffer_size,fp_out_groups);
          if(n_write!=buffer_size)
             SID_trap_error("buffer not written properly",ERROR_IO_WRITE);
       }
    }

    // Write group offsets
    int n_groups_written=0;
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          memcpy(buffer,offset_groups,sizeof(size_t)*n_groups_local);
          if(n_byte_offsets_groups==sizeof(unsigned int)){
             int i_buffer;
             for(i_buffer=0;i_buffer<n_groups_local;i_buffer++)
                buffer_int[i_buffer]=(int)buffer_size_t[i_buffer]; 
          }
          else if(n_byte_offsets_groups!=sizeof(int64_t))
             SID_trap_error("Invalid group offset byte-size (%d)",ERROR_LOGIC,n_byte_offsets_groups);
          buffer_size=n_byte_offsets_groups*n_groups_local;
       }
       SID_Bcast(&buffer_size,sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,      buffer_size,i_rank,SID.COMM_WORLD);
       if(SID.I_am_Master){
          int n_write;
          n_write=fwrite(buffer,1,buffer_size,fp_out_groups);
          if(n_write!=buffer_size)
             SID_trap_error("buffer not written properly",ERROR_IO_WRITE);
          n_groups_written+=n_write/n_byte_offsets_groups;
       }
    }

    // Write group subgroup counts
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          buffer_size=sizeof(int)*n_groups_local;
          memcpy(buffer,n_subgroups_group,buffer_size);
       }
       SID_Bcast(&buffer_size,sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,      buffer_size,i_rank,SID.COMM_WORLD);
       if(SID.I_am_Master){
          int n_write;
          n_write=fwrite(buffer,1,buffer_size,fp_out_groups);
          if(n_write!=buffer_size)
             SID_trap_error("buffer not written properly",ERROR_IO_WRITE);
       }
    }

    // Write subgroup sizes
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          buffer_size=sizeof(int)*n_subgroups_local;
          memcpy(buffer,size_subgroups,buffer_size);
       }
       SID_Bcast(&buffer_size,sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,      buffer_size,i_rank,SID.COMM_WORLD);
       if(SID.I_am_Master){
          int n_write;
          n_write=fwrite(buffer,1,buffer_size,fp_out_subgroups);
          if(n_write!=buffer_size)
             SID_trap_error("buffer not written properly",ERROR_IO_WRITE);
       }
    }

    // Write subgroup offsets
    int n_subgroups_written=0;
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          memcpy(buffer,offset_subgroups,sizeof(size_t)*n_subgroups_local);
          if(n_byte_offsets_subgroups==sizeof(unsigned int)){
             int i_buffer;
             for(i_buffer=0;i_buffer<n_subgroups_local;i_buffer++)
                buffer_int[i_buffer]=(int)buffer_size_t[i_buffer];
          }
          else if(n_byte_offsets_subgroups!=sizeof(int64_t))
             SID_trap_error("Invalid subgroup offset byte-size (%d)",ERROR_LOGIC,n_byte_offsets_subgroups);
          buffer_size=n_byte_offsets_subgroups*n_subgroups_local;
       }
       SID_Bcast(&buffer_size,sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,      buffer_size,i_rank,SID.COMM_WORLD);
       if(SID.I_am_Master){
          int n_write;
          n_write=fwrite(buffer,1,buffer_size,fp_out_subgroups);
          if(n_write!=buffer_size)
             SID_trap_error("buffer not written properly",ERROR_IO_WRITE);
          n_subgroups_written+=n_write/n_byte_offsets_subgroups;
       }
    }

    // Write particle IDs
    size_t n_particles_written=0;
    for(i_rank=0;i_rank<SID.n_proc;i_rank++){
       if(SID.My_rank==i_rank){
          memcpy(buffer,id_list_local,sizeof(size_t)*n_particles_local);
          if(n_byte_ids==sizeof(int)){
             int i_buffer;
             for(i_buffer=0;i_buffer<n_particles_local;i_buffer++)
                buffer_int[i_buffer]=(int)buffer_size_t[i_buffer];
          }
          else if(n_byte_ids!=sizeof(int64_t))
             SID_trap_error("Invalid ID byte-size (%d)",ERROR_LOGIC,n_byte_ids);
          buffer_size=n_byte_ids*n_particles_local;
       }
       SID_Bcast(&buffer_size,sizeof(int),i_rank,SID.COMM_WORLD);
       SID_Bcast(buffer,      buffer_size,i_rank,SID.COMM_WORLD);
       if(SID.I_am_Master){
          int n_write;
          n_write=fwrite(buffer,1,buffer_size,fp_out_particles);
          if(n_write!=buffer_size)
             SID_trap_error("buffer not written properly",ERROR_IO_WRITE);
          n_particles_written+=n_write/n_byte_ids;
       }
    }
    SID_log("No. of groups    written = %d",  SID_LOG_COMMENT,n_groups_written);
    SID_log("No. of subgroups written = %d",  SID_LOG_COMMENT,n_subgroups_written);
    SID_log("No. of particles written = %lld",SID_LOG_COMMENT,n_particles_written);
    SID_log("Done.",SID_LOG_CLOSE);
    if(SID.I_am_Master){
       fclose(fp_out_groups);
       fclose(fp_out_subgroups);
       fclose(fp_out_particles);
    }

    // Clean-up
    free_plist(&plist);
    SID_free(SID_FARG id_duplicates_local);
    SID_free(SID_FARG index_duplicates_local);
    SID_free(SID_FARG index_duplicates_local_index);
    SID_free(SID_FARG index_duplicates_group_local);
    SID_free(SID_FARG index_duplicates_subgroup_local);

    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);
  }

  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}

