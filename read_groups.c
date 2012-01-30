#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

void read_groups(char        *filename_groups_root,
                  int          i_file,
                  int          mode,
                  plist_info  *plist,
                  char        *catalog_name, ...){
   char    filename_cat[5];
   char    filename_groups[256];
   char    filename_subgroups[256];
   char    filename_ids[256];
   int     i,j,k,l;
   int     id_byte_size;
   int     n_ids_i;  
   int     i_rank;
   int     i_group;
   int     j_subgroup;
   int     i_subgroup;
   int     j_group;
   size_t  i_particle;
   size_t  j_particle;
   int     rank_offset_groups;
   int     rank_offset_subgroups;
   size_t  n_particles;
   int     n_particles_group;
   int     n_p_1,n_p_2;
   int    *n_groups_rank;
   int    *n_subgroups_rank;
   size_t *n_particles_rank;
   int    *global_group_index_local;
   int     n_groups;
   int     n_groups_i;
   int     n_subgroups;
   int     n_subgroups_local;
   int     n_subgroups_i;
   int    *n_subgroups_group;
   int    *read_index_subgroup;
   int     n_groups_in;
   int     n_subgroups_in;
   size_t  n_ids;
   int     n_groups_byte_length;
   int    *group_id;
   int    *group_length;
   int    *subgroup_length;
   int    *subgroup_offset;
   int    *group_offset;
   int     offset_0;
   size_t *input_id;
   int     id_in;
   long    record_length;
   char    parm_txt[256];
   double  expansion_factor;
   float   overdensity_vir;
   int     n_p_min;
   int     n_p_max;
   double *M_vir;
   int    *n_p_vir;
   float  *x_bound;
   float  *y_bound;
   float  *z_bound;
   double *x_COM;
   double *y_COM;
   double *z_COM;
   double *x_CODen;
   double *y_CODen;
   double *z_CODen;
   double *vx_COM;
   double *vy_COM;
   double *vz_COM;
   double *j_x;
   double *j_y;
   double *j_z;
   double *lambda;
   double *c_15;
   int     flag_read_groups;
   int     flag_read_ids;
   int     flag_read_subgroups;
   int     flag_read_MBP_ids_only=FALSE;
   int     flag_PHK_distribute;
   PHK_t  *keys_boundary=NULL;
   int     n_keys_boundary=0;
   va_list vargs;
   va_start(vargs,catalog_name);
 
   SID_profile_start("read_groups",SID_PROFILE_DEFAULT);
   SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
 
   // Set filenames
   sprintf(filename_cat,      "%03d",                   i_file);
   sprintf(filename_groups,   "%s_%s.catalog_groups",   filename_groups_root,filename_cat);
   sprintf(filename_subgroups,"%s_%s.catalog_subgroups",filename_groups_root,filename_cat);
   sprintf(filename_ids,      "%s_%s.catalog_particles",filename_groups_root,filename_cat);
 
   // Create flags from mode
   flag_read_groups   =TRUE;
   flag_read_ids      =TRUE;
   flag_read_subgroups=FALSE;
   flag_PHK_distribute=FALSE;
   if((check_mode_for_flag(mode,READ_GROUPS_ALL)        ||
       check_mode_for_flag(mode,READ_GROUPS_SUBGROUPS)) &&
      !check_mode_for_flag(mode,READ_GROUPS_NOSUBGROUPS))
     flag_read_subgroups=TRUE;
   if((check_mode_for_flag(mode,READ_GROUPS_NOIDS)))
     flag_read_ids=FALSE;
   if((check_mode_for_flag(mode,READ_GROUPS_MBP_IDS_ONLY)))
     flag_read_MBP_ids_only=TRUE;
   if((check_mode_for_flag(mode,READ_GROUPS_PEANOHILBERT)))
     flag_PHK_distribute=TRUE;

   // Process optional parameters
   int PHK_min_local_in;
   int PHK_max_local_in;
   if(flag_PHK_distribute){
      PHK_min_local_in=va_arg(vargs,int);
      PHK_max_local_in=va_arg(vargs,int);
   }
 
   // Read catalog file
   SID_log("Reading group data...",SID_LOG_OPEN|SID_LOG_TIMER);

   // Read header information ...
   int   header_size_groups   =0;
   int   header_size_subgroups=0;
   int   header_size_ids      =0;
   FILE *fp_groups   =NULL;
   FILE *fp_subgroups=NULL;
   FILE *fp_ids      =NULL;
   // ... from groups file...
   if((fp_groups=fopen(filename_groups,"r"))!=NULL){
      header_size_groups=sizeof(int);
      if(SID.I_am_Master)
         fread(&n_groups,sizeof(int),1,fp_groups);
      fclose(fp_groups);
      SID_Bcast(&n_groups,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
   }
   else
      SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_groups);

   // ... from subgroups file...
   if((fp_subgroups=fopen(filename_subgroups,"r"))!=NULL){
      header_size_subgroups=sizeof(int);
      if(SID.I_am_Master)
         fread(&n_subgroups,sizeof(int),1,fp_subgroups);
      fclose(fp_subgroups);
      SID_Bcast(&n_subgroups,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
   }
   else
      SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_subgroups);

   // ... from ids file...
   if((fp_ids=fopen(filename_ids,"r"))!=NULL){
      if(SID.I_am_Master){
         fread(&id_byte_size,sizeof(int),1,fp_ids);
         if(id_byte_size==sizeof(int)){
            fread(&n_ids_i,sizeof(int),1,fp_ids);
            n_ids=(size_t)n_ids_i;
         }
         else
            fread(&n_ids,sizeof(size_t),1,fp_ids);
      }
      fclose(fp_ids);
      SID_Bcast(&id_byte_size,sizeof(int),   MASTER_RANK,SID.COMM_WORLD);
      SID_Bcast(&n_ids,       sizeof(size_t),MASTER_RANK,SID.COMM_WORLD);
      header_size_ids=sizeof(int)+id_byte_size;
      if(id_byte_size>sizeof(size_t))
         SID_trap_error("Internal type size not sufficent for the ids in {%s}.",ERROR_LOGIC,filename_ids);
   }
   else
      SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_ids);

   // Store some header information
   ADaPS_store(&(plist->data),&n_groups,   "n_groups_all_%s",   ADaPS_SCALAR_INT,catalog_name);
   ADaPS_store(&(plist->data),&n_subgroups,"n_subgroups_all_%s",ADaPS_SCALAR_INT,catalog_name);
   if(flag_read_MBP_ids_only){
      size_t n_particiles_temp;
      n_particiles_temp=(size_t)n_groups;
      ADaPS_store(&(plist->data),&n_particiles_temp,"n_particles_all_%s",ADaPS_SCALAR_SIZE_T,catalog_name);
   }
   else
      ADaPS_store(&(plist->data),&n_particles,"n_particles_all_%s",ADaPS_SCALAR_SIZE_T,catalog_name);

   // Initialize domain decomposition ...
   char    filename_in_PHKs[MAX_FILENAME_LENGTH];
   int     n_bits_PHK;
   int    *PHK_groups;
   int    *PHK_subgroups;
   FILE   *fp_PHKs                =NULL;
   int     PHK_min_local          =0;
   int     PHK_max_local          =0;
   size_t  n_particles_local      =0;
   int     n_groups_local         =0;
   int    *read_index_group       =NULL;
   int    *read_seek_subgroup     =NULL;
   int    *read_seek_particles    =NULL;
   int    *storage_index_group    =NULL;
   int    *storage_index_subgroup =NULL;
   int    *storage_index_particles=NULL;
   int     i_buffer;
   int     n_buffer;
   int     n_buffer_max=1000;
   int     i_boundary;
   int     i_interior;
   int     n_groups_boundary;
   int     n_subgroups_boundary;
   int     n_groups_interior;
   n_groups_local   =0;
   n_subgroups_local=0;
   n_particles_local=0;

   // ... use precomputed Peano-Hilbert keys ...
   if(flag_PHK_distribute){
     SID_log("Reading PHKs...",SID_LOG_OPEN);

     // Open the file and check for success.  The first pointer will read the 
     //   PH key array while the second pointer will read the sorted particle counts.
     sprintf(filename_in_PHKs,"%s_%03d.catalog_PHKs",filename_groups_root,i_file);
     fp_PHKs=fopen(filename_in_PHKs,"r");
     if(fp_PHKs!=NULL){

        // Read the header values
        if(SID.I_am_Master){
           fread(&n_groups,   sizeof(int),   1,fp_PHKs);
           fread(&n_bits_PHK, sizeof(int),   1,fp_PHKs);
           fread(&n_particles,sizeof(size_t),1,fp_PHKs);
        }
        else
          fseeko(fp_PHKs,(off_t)(2*sizeof(int)+sizeof(size_t)),SEEK_SET);
        SID_Bcast(&n_groups,   sizeof(int),   MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(&n_bits_PHK, sizeof(int),   MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(&n_particles,sizeof(size_t),MASTER_RANK,SID.COMM_WORLD);
        SID_log("(%d groups, %lld particles; %d-bit keys)...",SID_LOG_CONTINUE,n_groups,n_particles,n_bits_PHK);

        // Assign key ranges to ranks on the fly as we read the file (buffer reads)
        //   The file contains triplets for each group of the sorted PHKs, PHK sort index and
        //   cumulative PHK-ranked particle count for the catalog.
        char   *buffer;
        int     unit_size;
        int     PHK_last;
        int     PHK_i;
        size_t  PHK_index_i;
        size_t  n_particles_i;
        int     n_groups_left;
        size_t  n_particles_base;
        size_t  n_particles_left;

        // Create read buffer
        unit_size=2*sizeof(int)+sizeof(size_t);
        buffer   =(char *)SID_malloc(n_buffer_max*unit_size);

        // Determine local key ranges if they were not passed to us ...
        if(PHK_min_local_in<0 || PHK_max_local_in<0){
           for(i_rank=0,i_group=0,i_buffer=n_buffer_max,n_groups_left=n_groups,n_particles_left=n_particles,n_particles_base=0,PHK_last=-1;i_rank<SID.n_proc && i_group<n_groups;i_rank++){

              // Each rank gets at least one halo (unless there are none left)
              int n_seek;
              if(i_rank==SID.My_rank){
                 // Buffer the read
                 n_seek=0;
                 if(i_buffer>=n_buffer_max){
                    n_buffer=MIN(n_buffer_max,n_groups-i_group);
                    n_seek+=n_buffer;
                    i_buffer=0;
                    fread(buffer,unit_size,n_buffer,fp_PHKs);
                 }
                 PHK_i        =((int    *)(&(buffer[i_buffer*unit_size]              )))[0];
                 PHK_index_i  =((int    *)(&(buffer[i_buffer*unit_size+1*sizeof(int)])))[0];
                 n_particles_i=((size_t *)(&(buffer[i_buffer*unit_size+2*sizeof(int)])))[0];
                 PHK_min_local    =PHK_last+1;
                 PHK_max_local    =PHK_i;
                 PHK_last         =PHK_i;
                 n_particles_local=n_particles_i-n_particles_base;
                 n_groups_local   =1;
                 i_group++;
                 i_buffer++;
                 while(i_group<n_groups){
                    // Buffer the read
                    if(i_buffer>=n_buffer_max){
                       n_buffer=MIN(n_buffer_max,n_groups-i_group);
                       n_seek+=n_buffer;
                       i_buffer=0;
                       fread(buffer,unit_size,n_buffer,fp_PHKs);
                    }
                    PHK_i        =((int    *)(&(buffer[i_buffer*unit_size]              )))[0];
                    PHK_index_i  =((int    *)(&(buffer[i_buffer*unit_size+1*sizeof(int)])))[0];
                    n_particles_i=((size_t *)(&(buffer[i_buffer*unit_size+2*sizeof(int)])))[0];
                    if(n_particles_local<(size_t)((float)(n_particles_left)/(float)(SID.n_proc-i_rank)) || PHK_i==PHK_max_local){
                       PHK_max_local    =PHK_i;
                       PHK_last         =PHK_i;
                       n_particles_local=n_particles_i-n_particles_base;
                       n_groups_local++;
                       i_group++;
                       i_buffer++;
                    }
                    else
                       break;
                 }
                 n_groups_left   -=n_groups_local;
                 n_particles_base+=n_particles_local;
                 n_particles_left-=n_particles_local;
              }

              // Bring the file pointers of all ranks up to date
              SID_Bcast(&n_seek,  sizeof(int),i_rank,SID.COMM_WORLD);
              SID_Bcast(&n_buffer,sizeof(int),i_rank,SID.COMM_WORLD);
              if(i_rank!=SID.My_rank && n_seek>0)
                 fseeko(fp_PHKs,(off_t)(n_seek*unit_size),SEEK_CUR);
              SID_Bcast(buffer,n_buffer*unit_size,i_rank,SID.COMM_WORLD);

              // Bring the particle and group counters up to date
              SID_Bcast(&n_particles_base,sizeof(size_t),i_rank,SID.COMM_WORLD);
              SID_Bcast(&n_groups_left,   sizeof(int),   i_rank,SID.COMM_WORLD);
              SID_Bcast(&n_particles_left,sizeof(size_t),i_rank,SID.COMM_WORLD);
              SID_Bcast(&i_group,         sizeof(int),   i_rank,SID.COMM_WORLD);
              SID_Bcast(&i_buffer,        sizeof(int),   i_rank,SID.COMM_WORLD);
              SID_Bcast(&PHK_last,        sizeof(int),   i_rank,SID.COMM_WORLD);
           }

           // If there are left-over ranks, they receieve nothing
           for(;i_rank<SID.n_proc;i_rank++){
              if(i_rank==SID.My_rank){
                 PHK_min_local    =PHK_DIM_SIZE(n_bits_PHK)*PHK_DIM_SIZE(n_bits_PHK)*PHK_DIM_SIZE(n_bits_PHK)+i_rank;
                 PHK_max_local    =PHK_min_local;
                 n_particles_local=0;
                 n_groups_local   =0;
              }
           }
        }
        // ... else, count the local number of particles and groups for the given range ...
        else{
           int PHK_last_rank;
           PHK_min_local =PHK_min_local_in;
           PHK_max_local =PHK_max_local_in;
           n_groups_local   =0;
           n_particles_local=0;

           // Check the given PHK ranges
           for(i_rank=0,PHK_last_rank=-1;i_rank<SID.n_proc;i_rank++){
              if(i_rank==SID.My_rank){
                 if(PHK_min_local!=(PHK_last_rank+1))
                    SID_trap_error("Keys not continuous at start of rank %d.",ERROR_LOGIC,SID.My_rank);
                 PHK_last_rank=PHK_max_local;
              }
              SID_Bcast(&PHK_last_rank,sizeof(int),i_rank,SID.COMM_WORLD);
           }

           for(i_group=0,i_buffer=n_buffer_max,n_particles_base=0;i_group<n_groups;i_group++,i_buffer++){
              // Buffer the read
              if(i_buffer>=n_buffer_max){
                 n_buffer=MIN(n_buffer_max,n_groups-i_group);
                 i_buffer=0;
                 if(SID.I_am_Master)
                    fread(buffer,unit_size,n_buffer,fp_PHKs);
                 SID_Bcast(buffer,n_buffer*unit_size,MASTER_RANK,SID.COMM_WORLD);
              }
              PHK_i        =((int    *)(&(buffer[i_buffer*unit_size]              )))[0];
              PHK_index_i  =((int    *)(&(buffer[i_buffer*unit_size+1*sizeof(int)])))[0];
              n_particles_i=((size_t *)(&(buffer[i_buffer*unit_size+2*sizeof(int)])))[0];
              if(PHK_i>=PHK_min_local && PHK_i<=PHK_max_local){
                 n_particles_local=n_particles_i-n_particles_base;
                 n_groups_local++;
              }
              else
                n_particles_base=n_particles_i;
           }

           // Determine if any groups or particles are unassigned
           int    n_groups_check;
           size_t n_particles_check;
           SID_Allreduce(&n_groups_local,   &n_groups_check,   1,SID_INT,   SID_SUM,SID.COMM_WORLD);
           SID_Allreduce(&n_particles_local,&n_particles_check,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
           n_particles_left=n_particles-n_particles_check;
           n_groups_left   =n_groups   -n_groups_check;
        } 

        // Sanity checks
        if(n_particles_left!=0)
           SID_trap_error("Particle number mismatch (%zd unassigned) after setting PHK ranges.",ERROR_LOGIC,n_particles_left);
        if(n_groups_left!=0)
           SID_trap_error("Group number mismatch (%d unassigned) after setting PHK ranges.",ERROR_LOGIC,n_groups_left);

        // Report results
        if(SID.n_proc>1){
           int PHK_min_i;
           int PHK_max_i;
           for(i_rank=0;i_rank<SID.n_proc;i_rank++){
              PHK_min_i=PHK_min_local;
              PHK_max_i=PHK_max_local;
              SID_Bcast(&PHK_min_i,sizeof(int),i_rank,SID.COMM_WORLD);
              SID_Bcast(&PHK_max_i,sizeof(int),i_rank,SID.COMM_WORLD);
              SID_log("PHK range for rank %4d: %d->%d",SID_LOG_COMMENT,i_rank,PHK_min_i,PHK_max_i);
           }
        }

        // Re-read the file to create local PHK and read index arrays
        int *PHK_group_temp;
        int *read_index_group_temp;
        fseeko(fp_PHKs,(off_t)(2*sizeof(int)+sizeof(size_t)),SEEK_SET);
        PHK_group_temp       =(int *)SID_malloc(sizeof(int)*n_groups_local);
        read_index_group_temp=(int *)SID_malloc(sizeof(int)*n_groups_local);
        for(i_group=0,j_group=0,i_buffer=n_buffer_max;i_group<n_groups;i_group++,i_buffer++){
           // Buffer the read
           if(i_buffer>=n_buffer_max){
              n_buffer=MIN(n_buffer_max,n_groups-i_group);
              i_buffer=0;
              if(SID.I_am_Master)
                 fread(buffer,unit_size,n_buffer,fp_PHKs);
              SID_Bcast(buffer,n_buffer*unit_size,MASTER_RANK,SID.COMM_WORLD);
           }
           PHK_i        =((int    *)(&(buffer[i_buffer*unit_size]              )))[0];
           PHK_index_i  =((int    *)(&(buffer[i_buffer*unit_size+1*sizeof(int)])))[0];
           // Store the key and read index in the local arrays if the current group a member of this rank.
           if(PHK_i>=PHK_min_local && PHK_i<=PHK_max_local){
              read_index_group_temp[j_group]=PHK_index_i;
              PHK_group_temp[j_group]       =PHK_i;
              j_group++;
           }
        }
        fclose(fp_PHKs);
        SID_free(SID_FARG buffer);

        // Sanity check
        if(j_group!=n_groups_local)
          SID_trap_error("Group counts don't make sense (ie. %d!=%d) after re-reading the PHK keys.",ERROR_LOGIC,j_group,n_groups_local);

        // Sort read indices
        size_t *read_index_group_temp_index;
        merge_sort(read_index_group_temp,(size_t)n_groups_local,&read_index_group_temp_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);

        // Sort index and PHK arrays (these are read-ordered, not storage-ordered)
        PHK_groups=(int *)SID_malloc(sizeof(int)*n_groups_local);
        read_index_group=(int *)SID_malloc(sizeof(int)*n_groups_local);
        for(i_group=0;i_group<n_groups_local;i_group++){
           PHK_groups[i_group]    =PHK_group_temp[read_index_group_temp_index[i_group]];        // Read-ordered
           read_index_group[i_group]=read_index_group_temp[read_index_group_temp_index[i_group]]; // Final read indices
        }
        SID_free(SID_FARG PHK_group_temp);
        SID_free(SID_FARG read_index_group_temp);
        SID_free(SID_FARG read_index_group_temp_index);

        // Determine which keys lie on the outside boundary of the local domain
        compute_PHK_boundary_keys(n_bits_PHK,PHK_min_local,PHK_max_local,&n_keys_boundary,&keys_boundary);

        // Determine how many groups have boundary keys and how many are in the interior
        for(i_group=0,n_groups_boundary=0,n_groups_interior=0;i_group<n_groups_local;i_group++){
          PHK_t key_temp;
          key_temp=(PHK_t)PHK_groups[i_group];
          if(is_a_member(&key_temp,keys_boundary,n_keys_boundary,SID_PHK_T))
             n_groups_boundary++;
          else
             n_groups_interior++;
        }

        // Create temporary arrays for boundary and interior groups
        int *PHK_boundary;
        int *PHK_interior;
	int *local_index_group_boundary;
        int *local_index_group_interior;
        PHK_boundary              =(int *)SID_malloc(sizeof(int)*n_groups_boundary); // read ordered
        local_index_group_boundary=(int *)SID_malloc(sizeof(int)*n_groups_boundary); // read indices
        PHK_interior              =(int *)SID_malloc(sizeof(int)*n_groups_interior); // read ordered
        local_index_group_interior=(int *)SID_malloc(sizeof(int)*n_groups_interior); // read indices
        for(i_group=0,i_boundary=0,i_interior=0;i_group<n_groups_local;i_group++){
           PHK_t key_temp;
           key_temp=(PHK_t)PHK_groups[i_group];
           if(is_a_member(&key_temp,keys_boundary,n_keys_boundary,SID_PHK_T)){
              PHK_boundary[i_boundary]              =PHK_groups[i_group];
              local_index_group_boundary[i_boundary]=i_group;
              i_boundary++;
           }
           else{
              PHK_interior[i_interior]              =PHK_groups[i_group];
              local_index_group_interior[i_interior]=i_group;
              i_interior++;
           }
        }

        // Sort boundary and interior keys separately
        size_t *PHK_boundary_index;
        size_t *PHK_interior_index;
        merge_sort(PHK_boundary,(size_t)n_groups_boundary,&PHK_boundary_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
        merge_sort(PHK_interior,(size_t)n_groups_interior,&PHK_interior_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);

        // Set the final storage indices.  Place boundary keys first and interiors after.
        storage_index_group=(int *)SID_malloc(sizeof(int)*n_groups_local);
        for(i_group=0,j_group=0;i_group<n_groups_boundary;i_group++,j_group++){
           PHK_groups[j_group]                                                         =PHK_boundary[PHK_boundary_index[i_group]];
           storage_index_group[local_index_group_boundary[PHK_boundary_index[i_group]]]=j_group;
        }
        for(i_group=0;i_group<n_groups_interior;i_group++,j_group++){
           PHK_groups[j_group]                                                         =PHK_interior[PHK_interior_index[i_group]];
           storage_index_group[local_index_group_interior[PHK_interior_index[i_group]]]=j_group;
        }

        // Sanity check
        if(j_group!=n_groups_local)
          SID_trap_error("Group counts don't make sense (ie. %d!=%d) after determining storage indices.",ERROR_LOGIC,j_group,n_groups_local);

        // Clean-up
        SID_free(SID_FARG PHK_boundary);
        SID_free(SID_FARG PHK_interior);
        SID_free(SID_FARG PHK_boundary_index);
        SID_free(SID_FARG PHK_interior_index);
        SID_free(SID_FARG local_index_group_boundary);
        SID_free(SID_FARG local_index_group_interior);

        // Store results
        ADaPS_store(&(plist->data),(void *)(&n_bits_PHK),        "n_bits_PHK_%s",          ADaPS_SCALAR_INT,catalog_name);
        ADaPS_store(&(plist->data),(void *)(&PHK_min_local),     "PHK_min_local_%s",       ADaPS_SCALAR_INT,catalog_name);
        ADaPS_store(&(plist->data),(void *)(&PHK_max_local),     "PHK_max_local_%s",       ADaPS_SCALAR_INT,catalog_name);
        ADaPS_store(&(plist->data),(void *)(&n_keys_boundary),   "n_keys_boundary_%s",     ADaPS_SCALAR_INT,catalog_name);
        ADaPS_store(&(plist->data),(void *)(&n_groups_boundary), "n_groups_boundary_%s",   ADaPS_SCALAR_INT,catalog_name);
        ADaPS_store(&(plist->data),(void *)(keys_boundary),      "keys_boundary_%s",       ADaPS_DEFAULT,   catalog_name);
        ADaPS_store(&(plist->data),(void *)(PHK_groups),         "PHK_groups_%s",          ADaPS_DEFAULT,   catalog_name);

        SID_log("Done.",SID_LOG_CLOSE);
     }
     // ... else if PHK file not opened ...
     else
        SID_trap_error("File {%s} not found.",ERROR_IO_OPEN,filename_in_PHKs);
   }
   // ... assign to ranks in rank-ordered blocks (default).
   else{
     SID_log("Initializing...",SID_LOG_OPEN);
     if((fp_groups=fopen(filename_groups,"r"))!=NULL){
        int    n_seek;
        int    group_length_i;
        size_t n_particles_left;
        if(SID.I_am_Master){
          fread(&n_groups,sizeof(int),1,fp_groups);
          for(i_group=0,n_particles=0;i_group<n_groups;i_group++){
             fread(&group_length_i,sizeof(int),1,fp_groups);
             n_particles+=(size_t)group_length_i;
          }
        }
        SID_Bcast(&n_groups,   sizeof(int),   MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(&n_particles,sizeof(size_t),MASTER_RANK,SID.COMM_WORLD);
        n_groups_local   =0;
        n_particles_local=0;
        n_seek           =0;
        fseeko(fp_groups,(off_t)(sizeof(int)),SEEK_SET);
        for(i_rank=0,i_group=0,n_particles_left=n_particles;i_rank<SID.n_proc && n_particles_left>0;i_rank++){
           if(n_seek>0)
              fseeko(fp_groups,(off_t)(n_seek*sizeof(int)),SEEK_CUR);
           if(i_rank==SID.My_rank){
              n_seek=0;
              while(i_group<n_groups && n_particles_local<(size_t)((float)(n_particles_left)/(float)(SID.n_proc-i_rank))){
                 fread(&group_length_i,sizeof(int),1,fp_groups);
                 n_groups_local++;
                 n_particles_local+=(size_t)group_length_i;
                 n_seek++;
                 i_group++;
              }
              n_particles_left-=n_particles_local;
           }
           SID_Bcast(&n_particles_left,sizeof(size_t),i_rank,SID.COMM_WORLD);
           SID_Bcast(&n_seek,          sizeof(int),   i_rank,SID.COMM_WORLD);
           SID_Bcast(&i_group,         sizeof(int),   i_rank,SID.COMM_WORLD);
           if(i_rank==SID.My_rank)
              n_seek=0;
        }
        // If there are left-over ranks, they receieve nothing
        for(;i_rank<SID.n_proc;i_rank++){
           if(i_rank==SID.My_rank){
              n_particles_local=0;
              n_groups_local   =0;
           }
        }
        calc_sum_global(&n_particles_local,&n_particles,1,SID_SIZE_T,CALC_MODE_DEFAULT,SID.COMM_WORLD);
        fclose(fp_groups);
     }
     else{
        n_groups   =0;
        n_particles=0;
     }

     // Create group read and storage arrays
     read_index_group   =(int *)SID_malloc(sizeof(int)*n_groups_local);
     storage_index_group=(int *)SID_malloc(sizeof(int)*n_groups_local);
     for(i_rank=0,i_group=0;i_rank<SID.n_proc;i_rank++){
        if(i_rank==SID.My_rank){
           for(j_group=0;j_group<n_groups_local;i_group++,j_group++){
             read_index_group[j_group]   =i_group;
             storage_index_group[j_group]=j_group;
           }
        }
        SID_Bcast(&i_group,sizeof(int),i_rank,SID.COMM_WORLD);
     }
     SID_log("(%d groups, %lld particles)...Done.",SID_LOG_CLOSE,n_groups,n_particles);
   }

   // This array will be needed in the future to connect groups in memory to their positions on disk
   //   It is created by placing the read-ordered read indices into storage-order
   int *read_index_group_store;
   read_index_group_store=(int *)SID_malloc(sizeof(int)*n_groups_local);
   for(i_group=0;i_group<n_groups_local;i_group++)
     read_index_group_store[storage_index_group[i_group]]=read_index_group[i_group];
   ADaPS_store(&(plist->data),(void *)(read_index_group_store),"file_index_groups_%s",ADaPS_DEFAULT,catalog_name);

   // Consistancy checks
   int    n_groups_check;
   size_t n_particles_check;
   SID_Allreduce(&n_groups_local,&n_groups_check,1,SID_INT,SID_SUM,SID.COMM_WORLD);
   if(n_groups_check!=n_groups)
      SID_trap_error("Inconsistant group count (%d!=%d) after initializing domain decomposition.",ERROR_LOGIC,n_groups_check,n_groups);
   SID_Allreduce(&n_particles_local,&n_particles_check,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
   if(n_particles_check!=n_particles)
      SID_trap_error("Inconsistant particle count (%zd!=%zd) after initializing domain decomposition.",ERROR_LOGIC,n_particles_check,n_particles);
 
   // Read the files ...
   if(flag_read_groups||flag_read_subgroups){
 
      // Read the group file ...
      if(flag_read_groups){
         SID_log("Reading group file {%s}...",SID_LOG_OPEN,filename_groups);
         if((fp_groups=fopen(filename_groups,"r"))!=NULL)
            fseeko(fp_groups,(off_t)header_size_groups,SEEK_SET);
         else
           SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_groups);
         if(n_groups>0){
            // ... allocate arrays ...
            int *buffer;
            group_length     =(int *)SID_malloc(sizeof(int)*n_groups_local);
            group_offset     =(int *)SID_malloc(sizeof(int)*n_groups_local);
            n_subgroups_group=(int *)SID_malloc(sizeof(int)*n_groups_local);

            // These variables are used to check for indexing inconsistancies
            int  i_test,j_test;
            int *test;
            test=(int *)SID_malloc(sizeof(int)*n_groups_local);

            // ... read local group lengths ...
            buffer=(int *)SID_malloc(sizeof(int)*n_buffer_max);
            for(i_test=0;i_test<n_groups_local;i_test++) test[i_test]=0;
            for(i_group=0,j_group=0,i_buffer=n_buffer_max;i_group<n_groups;i_group++,i_buffer++){
               // Buffer the read
               if(i_buffer>=n_buffer_max){
                  i_buffer=0;
                  n_buffer=MIN(n_buffer_max,n_groups-i_group);
                  if(SID.I_am_Master)
                     fread(buffer,sizeof(int),n_buffer,fp_groups);
                  SID_Bcast(buffer,n_buffer*sizeof(int),MASTER_RANK,SID.COMM_WORLD);
               }
               if(j_group<n_groups_local){
                  if(read_index_group[j_group]==i_group){
                     group_length[storage_index_group[j_group]]=buffer[i_buffer];
                     test[storage_index_group[j_group]]++;
                     j_group++;
                  }
               }
            }
            if(j_group!=n_groups_local)
              SID_trap_error("Group counts don't make sense (ie. %d!=%d) after reading group lengths.",ERROR_LOGIC,j_group,n_groups_local);

            // Check for indexing inconsistancies
            for(i_group=0;i_group<n_groups_local;i_group++){
               if(test[i_group]!=1)
                  SID_trap_error("Group length indexing error: %5d %5d\n",i_group,test[i_group]);
            }

            // Count the number of particles in the boundary keys
            if(flag_PHK_distribute){
               size_t n_particles_boundary;
               for(i_group=0,n_particles_boundary=0;i_group<n_groups_boundary;i_group++)
                  n_particles_boundary+=group_length[i_group];
               ADaPS_store(&(plist->data),(void *)(&n_particles_boundary),"n_particles_boundary_%s",ADaPS_SCALAR_SIZE_T,catalog_name);
            }

            // ... read group offsets ...
            for(i_test=0;i_test<n_groups_local;i_test++) test[i_test]=0;
            for(i_group=0,j_group=0,i_buffer=n_buffer_max;i_group<n_groups;i_group++,i_buffer++){
               // Buffer the read
               if(i_buffer>=n_buffer_max){
                  i_buffer=0;
                  n_buffer=MIN(n_buffer_max,n_groups-i_group);
                  if(SID.I_am_Master)
                     fread(buffer,sizeof(int),n_buffer,fp_groups);
                  SID_Bcast(buffer,n_buffer*sizeof(int),MASTER_RANK,SID.COMM_WORLD);
               }
               if(j_group<n_groups_local){
                  if(read_index_group[j_group]==i_group){
                     group_offset[storage_index_group[j_group]]=buffer[i_buffer];
                     test[storage_index_group[j_group]]++;
                     j_group++;
                  }
               }
            }
            if(j_group!=n_groups_local)
              SID_trap_error("Group counts don't make sense (ie. %d!=%d) after reading group offsets.",ERROR_LOGIC,j_group,n_groups_local);

            // Check for indexing inconsistancies
            for(i_group=0;i_group<n_groups_local;i_group++){
               if(test[i_group]!=1)
                  SID_trap_error("Group offset indexing error: %5d %5d\n",i_group,test[i_group]);
            }

            // ... read number of subgroups per group ...
            int *read_index_subgroup_temp;
            read_seek_subgroup      =(int *)SID_calloc(sizeof(int)*(n_groups_local+1)); // Must be initialized to zero; one extra needed for end-block syncing
            read_index_subgroup_temp=(int *)SID_malloc(sizeof(int)*n_groups_local);
            for(i_test=0;i_test<n_groups_local;i_test++) test[i_test]=0;
            for(i_group=0,j_group=0,i_buffer=n_buffer_max,i_subgroup=0,n_subgroups_local=0;i_group<n_groups;i_group++,i_buffer++){
               // Buffer the read
               if(i_buffer>=n_buffer_max){
                  i_buffer=0;
                  n_buffer=MIN(n_buffer_max,n_groups-i_group);
                  if(SID.I_am_Master)
                     fread(buffer,sizeof(int),n_buffer,fp_groups);
                  SID_Bcast(buffer,n_buffer*sizeof(int),MASTER_RANK,SID.COMM_WORLD);
               }
               if(j_group<n_groups_local){
                  if(read_index_group[j_group]==i_group){
                     n_subgroups_local                                     +=buffer[i_buffer];
                     n_subgroups_group[storage_index_group[j_group]]        =buffer[i_buffer];
                     read_index_subgroup_temp[storage_index_group[j_group]] =i_subgroup;
                     test[storage_index_group[j_group]]++;
                     j_group++;
                  }
                  else 
                     read_seek_subgroup[j_group]+=buffer[i_buffer];
               }
               else
                  read_seek_subgroup[j_group]+=buffer[i_buffer];
               i_subgroup+=buffer[i_buffer];
            }
            
            SID_free(SID_FARG buffer);
            if(j_group!=n_groups_local)
              SID_trap_error("Group counts don't make sense (ie. %d!=%d) after reading no. of subgroups per group.",ERROR_LOGIC,j_group,n_groups_local);
            if(i_subgroup!=n_subgroups)
              SID_trap_error("Subgroup counts don't make sense (ie. %d!=%d) after reading no. of subgroups per group.",ERROR_LOGIC,i_subgroup,n_subgroups);

            // Check for indexing inconsistancies
            for(i_group=0;i_group<n_groups_local;i_group++){
               if(test[i_group]!=1)
                  SID_trap_error("Group offset indexing error: %5d %5d\n",i_group,test[i_group]);
            }
            SID_free(SID_FARG test);

            // Subgroup count consistancy check
            int n_subgroups_check;
            calc_sum(        n_subgroups_group,&n_subgroups_local,n_groups_local,SID_INT,CALC_MODE_DEFAULT);
            calc_sum_global(&n_subgroups_local,&n_subgroups_check,1,             SID_INT,CALC_MODE_DEFAULT,SID.COMM_WORLD);
            if(n_subgroups_check!=n_subgroups)
               SID_trap_error("Inconsistant subgroup count (%d!=%d) after reading group file.",ERROR_LOGIC,n_subgroups_check,n_subgroups);
            if(i_subgroup!=n_subgroups)
              SID_trap_error("Local subgroup counts don't make sense (ie. %d!=%d) after reading no. of subgroups per group.",ERROR_LOGIC,i_subgroup,n_subgroups_local);

            // Finish computing subgroup read indices
            read_index_subgroup=(int *)SID_malloc(sizeof(int)*n_subgroups_local);
            for(i_group=0,i_subgroup=0;i_group<n_groups_local;i_group++){
               if(n_subgroups_group[i_group]>0){ // Some groups have no subgroups.  We don't want to increament i_subgroup in that case.
                 for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];i_subgroup++,j_subgroup++)
                    read_index_subgroup[i_subgroup]=read_index_subgroup_temp[i_group]+j_subgroup;
               }
            }
            SID_free(SID_FARG read_index_subgroup_temp);

            // Create subgroup and particle storage indices ...
            int *temp_array;
            temp_array             =(int *)SID_malloc(sizeof(int)*n_groups_local);
            storage_index_subgroup =(int *)SID_malloc(sizeof(int)*n_groups_local);
            storage_index_particles=(int *)SID_malloc(sizeof(int)*n_groups_local);

            // ... subgroups ...
            temp_array[0]=0;
            for(i_group=1;i_group<n_groups_local;i_group++)
               temp_array[i_group]=temp_array[i_group-1]+n_subgroups_group[i_group-1];
            for(i_group=0;i_group<n_groups_local;i_group++)
               storage_index_subgroup[i_group]=temp_array[storage_index_group[i_group]];

            // ... particles ...
            for(i_group=1;i_group<n_groups_local;i_group++)
               temp_array[i_group]=temp_array[i_group-1]+group_length[i_group-1];
            for(i_group=0;i_group<n_groups_local;i_group++)
               storage_index_particles[i_group]=temp_array[storage_index_group[i_group]];
            SID_free(SID_FARG temp_array);

            // Count the number of subgroups in the boundary region if we are using PHK decomposition
            if(flag_PHK_distribute){
               n_subgroups_boundary=0;
               for(i_group=0;i_group<n_groups_boundary;i_group++)
                  n_subgroups_boundary+=n_subgroups_group[i_group];
               ADaPS_store(&(plist->data),(void *)(&n_subgroups_boundary),"n_subgroups_boundary_%s",ADaPS_SCALAR_INT,catalog_name);
            }

            // Create subgroup PHK indices
            if(flag_PHK_distribute){
               PHK_subgroups=(int *)SID_malloc(sizeof(int)*n_subgroups_local);
               for(i_group=0,i_subgroup=0;i_group<n_groups_local;i_group++){
                  if(n_subgroups_group[i_group]>0){ // Some groups have no subgroups.  We don't want to increament i_subgroup in that case.
                    for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];i_subgroup++,j_subgroup++)
                       PHK_subgroups[i_subgroup]=PHK_groups[i_group];
                  }
               }
               ADaPS_store(&(plist->data),PHK_subgroups,"PHK_subgroups_%s",ADaPS_DEFAULT,catalog_name);
            }


            // Report demcomposition results
            if(SID.n_proc>1){
               SID_log("Results of domain decomposition:",SID_LOG_OPEN);
               int    n_groups_report;
               int    n_subgroups_report;
               size_t n_particles_report;
               for(i_rank=0;i_rank<SID.n_proc;i_rank++){
                  n_groups_report   =n_groups_local;
                  n_subgroups_report=n_subgroups_local;
                  n_particles_report=n_particles_local;
                  SID_Bcast(&n_groups_report,   sizeof(int),   i_rank,SID.COMM_WORLD);
                  SID_Bcast(&n_subgroups_report,sizeof(int),   i_rank,SID.COMM_WORLD);
                  SID_Bcast(&n_particles_report,sizeof(size_t),i_rank,SID.COMM_WORLD);
                  SID_log("(i_rank,n_groups,n_subgroups,n_particles)=(%3d,%7d,%7d,%9lld)",SID_LOG_COMMENT,
                          i_rank,n_groups_report,n_subgroups_report,n_particles_report);
               }
               SID_log("Done.",SID_LOG_CLOSE);
            }

            // Store results
            ADaPS_store(&(plist->data),group_length,       "n_particles_group_%s",    ADaPS_DEFAULT,   catalog_name);
            ADaPS_store(&(plist->data),group_offset,       "particle_offset_group_%s",ADaPS_DEFAULT,   catalog_name);
            ADaPS_store(&(plist->data),n_subgroups_group,  "n_subgroups_group_%s",    ADaPS_DEFAULT,   catalog_name);
            ADaPS_store(&(plist->data),read_index_subgroup,"file_index_subgroups_%s", ADaPS_DEFAULT,   catalog_name);
    
            SID_log("Done. (%d groups)",SID_LOG_CLOSE,n_groups);
         }
         else
            SID_log("NO GROUPS TO READ!",SID_LOG_CLOSE);
         fclose(fp_groups);
      }
 
      // Store some decomposition header information
      ADaPS_store(&(plist->data),&n_groups_local,   "n_groups_%s",   ADaPS_SCALAR_INT,catalog_name);
      ADaPS_store(&(plist->data),&n_subgroups_local,"n_subgroups_%s",ADaPS_SCALAR_INT,catalog_name);
      if(flag_read_MBP_ids_only){
         size_t n_particiles_temp;
         n_particiles_temp=(size_t)n_groups_local;
         ADaPS_store(&(plist->data),&n_particiles_temp,"n_particles_%s",ADaPS_SCALAR_SIZE_T,catalog_name);
      }
      else
         ADaPS_store(&(plist->data),(void *)(&(n_particles_local)),"n_particles_%s",    ADaPS_SCALAR_SIZE_T,catalog_name);
 
      // Read the subgroups file...
      if(flag_read_subgroups){
         SID_log("Reading subgroup file {%s}...",SID_LOG_OPEN,filename_subgroups);
         if((fp_subgroups=fopen(filename_subgroups,"r"))!=NULL)
            fseeko(fp_subgroups,(off_t)header_size_subgroups,SEEK_SET);
         else
           SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_subgroups);

         if(n_subgroups>0){
            int n_subgroups_check;
            int i_subgroup;
            int j_subgroup;
            int max_group_size;
 
            // Allocate arrays
            subgroup_length=(int *)SID_calloc(sizeof(int)*n_subgroups_local);
            subgroup_offset=(int *)SID_malloc(sizeof(int)*n_subgroups_local);

            // These variables are used to check for indexing inconsistancies
            int  i_test,j_test;
            int *test;
            test=(int *)SID_malloc(sizeof(int)*n_subgroups_local);

            // ... read local subgroup lengths ...
            int index_group;
            int index_subgroup;
            for(i_test=0;i_test<n_subgroups_local;i_test++) test[i_test]=0;
            for(i_group=0,j_group=0;i_group<n_groups && j_group<n_groups_local;i_group++){
               if(read_index_group[j_group]==i_group){
                  index_group   =storage_index_group[j_group];
                  index_subgroup=storage_index_subgroup[j_group];
                  fseeko(fp_subgroups,(off_t)(read_seek_subgroup[j_group]*sizeof(int)),SEEK_CUR);
                  fread(&(subgroup_length[index_subgroup]),sizeof(int),n_subgroups_group[index_group],fp_subgroups);
                  for(i_test=index_subgroup,j_test=0;j_test<n_subgroups_group[index_group];i_test++,j_test++) 
                     test[i_test]++;
                  j_group++;
               }
            }

            // Check for indexing inconsistancies
            for(i_group=0,i_subgroup=0;i_group<n_groups_local;i_group++){
               for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];i_subgroup++,j_subgroup++){
                  if(test[i_subgroup]!=1)
                     SID_trap_error("Subgroup length indexing error: %5d %5d %5d\n",i_group,j_subgroup,test[i_subgroup]);
               }
            }

            // ... (make sure all ranks are synced to the end of the subgroup lengths) ...
            fseeko(fp_subgroups,(off_t)(read_seek_subgroup[j_group]*sizeof(int)),SEEK_CUR);
            if(j_group!=n_groups_local)
              SID_trap_error("Group counts don't make sense (ie. %d!=%d) after reading group lengths.",ERROR_LOGIC,j_group,n_groups_local);

            // ... read local subgroup offsets ...
            for(i_test=0;i_test<n_subgroups_local;i_test++) test[i_test]=0;
            for(i_group=0,j_group=0;i_group<n_groups && j_group<n_groups_local;i_group++){
               if(read_index_group[j_group]==i_group){
                  index_group   =storage_index_group[j_group];
                  index_subgroup=storage_index_subgroup[j_group];
                  fseeko(fp_subgroups,(off_t)(read_seek_subgroup[j_group]*sizeof(int)),SEEK_CUR);
                  fread(&(subgroup_offset[index_subgroup]),sizeof(int),n_subgroups_group[index_group],fp_subgroups);
                  for(i_test=index_subgroup,j_test=0;j_test<n_subgroups_group[index_group];i_test++,j_test++) 
                     test[i_test]++;
                  // Sanity Check
                  for(i_subgroup=0;i_subgroup<n_subgroups_group[index_group];i_subgroup++){
                     if((subgroup_offset[index_subgroup+i_subgroup]+subgroup_length[index_subgroup+i_subgroup])>(group_offset[index_group]+group_length[index_group])){
                        SID_trap_error("Subgroup {%d}'s ID list over-runs group {%d}'s ID list.",ERROR_LOGIC,index_subgroup+i_subgroup,index_group);
                     }
                  }
                  j_group++;
               }
            }
            // ... (make sure all ranks are synced to the end of the subgroup lengths) ...
            fseeko(fp_subgroups,(off_t)(read_seek_subgroup[j_group]*sizeof(int)),SEEK_CUR);
            if(j_group!=n_groups_local)
              SID_trap_error("Group counts don't make sense (ie. %d!=%d) after reading group lengths.",ERROR_LOGIC,j_group,n_groups_local);

            // Check for indexing inconsistancies
            for(i_group=0,i_subgroup=0;i_group<n_groups_local;i_group++){
               for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];i_subgroup++,j_subgroup++){
                  if(test[i_subgroup]!=1)
                     SID_trap_error("Subgroup length indexing error: %5d %5d %5d",i_group,j_subgroup,test[i_subgroup]);
               }
            }
            SID_free(SID_FARG test);

            // ... sanity check ...
            for(i_group=0,i_subgroup=0;i_group<n_groups_local && i_subgroup<n_subgroups_local;i_subgroup+=n_subgroups_group[i_group],i_group++){
               if(n_subgroups_group[i_group]>0){
                 if(group_offset[i_group]!=subgroup_offset[i_subgroup])
                    SID_trap_error("The first subgroup (%d) in a group (%d) does not share the group's particle offset (ie %d!=%d).",ERROR_LOGIC,
                                   i_subgroup,i_group,group_offset[i_group],subgroup_offset[i_subgroup]);
               }
            }

            // Store results
            ADaPS_store(&(plist->data),subgroup_length,"n_particles_subgroup_%s",    ADaPS_DEFAULT,   catalog_name);
            ADaPS_store(&(plist->data),subgroup_offset,"particle_offset_subgroup_%s",ADaPS_DEFAULT,   catalog_name);
            fclose(fp_subgroups);
            SID_log("Done. (%ld subgroups)",SID_LOG_CLOSE,n_subgroups);
         }
         else
            SID_log("NO SUBGROUPS TO READ!",SID_LOG_CLOSE);
      }

      // Read the particle IDs file ...
      if(flag_read_ids){
         SID_log("Reading IDs file {%s}...",SID_LOG_OPEN,filename_ids);
         if((fp_ids=fopen(filename_ids,"r"))!=NULL)
            fseeko(fp_ids,(off_t)header_size_ids,SEEK_SET);
         else
           SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_ids);
         if(n_ids>0){

            // Create the array that will hold the IDs
            SID_log("%d byte IDs...",SID_LOG_CONTINUE,id_byte_size);
            if(flag_read_MBP_ids_only)
               input_id=(size_t *)SID_malloc(sizeof(size_t)*n_groups_local);
            else
               input_id=(size_t *)SID_malloc(sizeof(size_t)*n_particles_local); //**

            // These variables are used to check for indexing inconsistancies
            int  i_test,j_test;
            int *test;
            test=(int *)SID_calloc(sizeof(int)*n_particles_local);

            // Create read buffer
            int     max_group_length;
            size_t *buffer;
            void   *buffer_in;
            calc_max(group_length,&max_group_length,n_groups_local,SID_INT,CALC_MODE_DEFAULT);
            buffer_in=SID_malloc(id_byte_size*max_group_length);
            if(id_byte_size==sizeof(int))
               buffer=(size_t *)SID_malloc(sizeof(size_t)*max_group_length);
            else
               buffer=(size_t *)buffer_in;

            // Read the IDs
            int   index_group;
            int   index_particles;
            int   index_last;
            for(i_test=0;i_test<n_particles_local;i_test++) test[i_test]=0;
            fseeko(fp_ids,(off_t)header_size_ids,SEEK_SET);
            for(i_group=0,j_group=0,index_last=0;i_group<n_groups && j_group<n_groups_local;i_group++){
               if(read_index_group[j_group]==i_group){
                  index_group    =storage_index_group[j_group];
                  index_particles=storage_index_particles[j_group];
                  fseeko(fp_ids,(off_t)((group_offset[index_group]-index_last)*id_byte_size),SEEK_CUR);
                  fread(buffer_in,id_byte_size,group_length[index_group],fp_ids);
                  index_last=group_offset[index_group]+group_length[index_group];
                  if(buffer!=buffer_in){
                     for(i_buffer=0;i_buffer<group_length[index_group];i_buffer++)
                       buffer[i_buffer]=(size_t)(((int *)buffer_in)[i_buffer]);
                  }
                  if(flag_read_MBP_ids_only){
                     input_id[index_group]=buffer[0];
                     test[index_group]++;
                  }
                  else{
                     memcpy(&(input_id[index_particles]),buffer,sizeof(size_t)*group_length[index_group]);
                     for(i_test=index_particles,j_test=0;j_test<group_length[index_group];i_test++,j_test++)
                        test[i_test]++;
                  }
                  j_group++;
               }
            }
            if(j_group!=n_groups_local)
              SID_trap_error("Group counts don't make sense (ie. %d!=%d) after reading group IDs.",ERROR_LOGIC,j_group,n_groups_local);

            // Check for indexing inconsistancies
            if(flag_read_MBP_ids_only){
               for(i_group=0,i_particle=0;i_group<n_groups_local;i_group++,i_particle++){
                  if(test[i_particle]!=1)
                     SID_trap_error("Particle ID indexing error: %5d %5d",i_group,test[i_particle]);
               }
            }
            else{
               for(i_group=0,i_particle=0;i_group<n_groups_local;i_group++){
                  for(j_particle=0;j_particle<group_length[i_group];i_particle++,j_particle++){ 
                     if(test[i_particle]!=1)
                        SID_trap_error("Particle ID indexing error: %5d %5d",i_group,test[i_particle]);
                  }
               }
               if(i_particle!=n_particles_local)
                  SID_trap_error("Particle counts don't make sense (ie. %d!=%d) after sanity check of IDs indexing.",ERROR_LOGIC,i_particle,n_particles_local);
            }
            SID_free(SID_FARG test);

            // Recompute local particle offsets so they correspond to what's in RAM
            //   (This array is already stored in the ADaPS structure.  Changes here will carry through.)
            // ... compute subgroup offsets relative to their groups ...
            for(i_group=0,i_subgroup=0;i_group<n_groups_local && i_subgroup<n_subgroups_local;i_group++){
               for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];j_subgroup++,i_subgroup++){
                  subgroup_offset[i_subgroup]-=group_offset[i_group];
               }
            }

            // ... groups ...
            group_offset[0]=0;
            for(i_group=1;i_group<n_groups_local;i_group++)
               group_offset[i_group]=group_offset[i_group-1]+group_length[i_group-1];

            // ... finalize subgroup offsets ...
            for(i_group=0,i_subgroup=0;i_group<n_groups_local && i_subgroup<n_subgroups_local;i_group++){
               for(j_subgroup=0;j_subgroup<n_subgroups_group[i_group];j_subgroup++,i_subgroup++){
                  subgroup_offset[i_subgroup]+=group_offset[i_group];
               }
            }

            // Clean-up
            if(buffer!=buffer_in)
               SID_free(SID_FARG buffer);
            SID_free(SID_FARG buffer_in);

            SID_log("Done. (%ld particles)",SID_LOG_CLOSE,n_particles);
         }
         else{
            input_id=NULL;
            SID_log("NO IDS TO READ!",SID_LOG_CLOSE);
         }

         // Store results
         ADaPS_store(&(plist->data),(void *)(input_id),"particle_ids_%s",ADaPS_DEFAULT,catalog_name);

         fclose(fp_ids);
      }
   }

   // Clean-up
   SID_free(SID_FARG read_index_group);
   SID_free(SID_FARG read_seek_subgroup);
   SID_free(SID_FARG read_seek_particles);
   SID_free(SID_FARG storage_index_group);
   SID_free(SID_FARG storage_index_subgroup);
   SID_free(SID_FARG storage_index_particles);
   va_end(vargs);
   
   SID_log("Done.",SID_LOG_CLOSE);
   SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1); 
   SID_profile_stop(SID_PROFILE_DEFAULT);
}
