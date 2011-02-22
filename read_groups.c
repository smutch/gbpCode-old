#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpHalos.h>

// To do: add group_count_byte_length
//        add big_int to SID.h

void read_groups(char        *filename_groups_root,
                 int          i_file,
                 int          mode,
                 plist_info  *plist,
                 char        *catalog_name){
  char    filename_cat[5];
  char    filename_groups[256];
  char    filename_subgroups[256];
  char    filename_ids[256];
  FILE   *fp;
  int     i,j,k,l;
  size_t  i_particle;
  int     id_byte_size;
  int     n_ids_i;  
  int     i_rank;
  int     i_group;
  int     j_group;
  int     rank_offset;
  size_t  n_particles;
  size_t  n_particles_left;
  int     n_particles_group;
  int     n_p_1,n_p_2;
  int    *n_groups_rank;
  int    *n_subgroups_rank;
  size_t *n_particles_rank;
  int     n_groups;
  int     n_groups_i;
  int     n_subgroups;
  int     n_subgroups_i;
  int    *n_subgroups_group;
  int     n_groups_in;
  int     n_subgroups_in;
  size_t  n_ids;
  int     n_groups_byte_length;
  int    *group_id;
  int    *group_length;
  int    *subgroup_length;
  int    *subgroup_offset;
  int    *group_offsets;
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

  int test_max;

  SID_profile_start("read_groups",SID_PROFILE_DEFAULT);

  // Set filenames
  sprintf(filename_cat,      "%03d",                   i_file);
  sprintf(filename_groups,   "%s_%s.catalog_groups",   filename_groups_root,filename_cat);
  sprintf(filename_subgroups,"%s_%s.catalog_subgroups",filename_groups_root,filename_cat);
  sprintf(filename_ids,      "%s_%s.catalog_particles",filename_groups_root,filename_cat);

  // Create flags from mode
  flag_read_groups   =TRUE;
  flag_read_ids      =TRUE;
  flag_read_subgroups=FALSE;
  if((check_mode_for_flag(mode,READ_GROUPS_ALL)           ||
      check_mode_for_flag(mode,READ_GROUPS_SUBGROUPS))    &&
     !check_mode_for_flag(mode,READ_GROUPS_NOSUBGROUPS))
    flag_read_subgroups=TRUE;
  if((check_mode_for_flag(mode,READ_GROUPS_NOIDS)))
    flag_read_ids=FALSE;

  // Read catalog file
  SID_log("Reading group data...",SID_LOG_OPEN|SID_LOG_TIMER);
  
  // group file
  if(flag_read_groups){
    SID_log("Reading group catalogs {%s}...",SID_LOG_OPEN,filename_groups);
    if((fp=fopen(filename_groups,"r"))!=NULL){
      if(SID.I_am_Master)
        fread(&n_groups,sizeof(int),1,fp);
      else
        fseeko(fp,sizeof(int),SEEK_CUR);
      SID_Bcast(&n_groups,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
      
      // Check that >0 groups are being loaded
      n_particles=0;
      n_subgroups=0;
      group_length     =(int    *)SID_malloc(sizeof(int   )*n_groups);
      group_offsets    =(int    *)SID_malloc(sizeof(int   )*n_groups);
      n_subgroups_group=(int    *)SID_malloc(sizeof(int   )*n_groups);
      n_groups_rank    =(int    *)SID_malloc(sizeof(int   )*SID.n_proc);
      n_subgroups_rank =(int    *)SID_malloc(sizeof(int   )*SID.n_proc);
      n_particles_rank =(size_t *)SID_malloc(sizeof(size_t)*SID.n_proc);
      for(i_rank=0;i_rank<SID.n_proc;i_rank++){
        n_groups_rank[i_rank]   =0;
        n_subgroups_rank[i_rank]=0;
        n_particles_rank[i_rank]=0;
      }
      if(n_groups>0){
        // Read catalog indices
        
        // Figure-out how to distribute everything to the ranks
        if(SID.I_am_Master){
          fread(group_length,      sizeof(int),n_groups,fp);
          fread(group_offsets,     sizeof(int),n_groups,fp);
          fread(n_subgroups_group, sizeof(int),n_groups,fp);
          // Determine how many groups go on each rank
          for(i_group=0,n_particles=0;i_group<n_groups;i_group++){
            n_subgroups+=n_subgroups_group[i_group];
            n_particles+=group_length[i_group];
          }
          for(i_rank=0,i_group=0,n_particles_left=n_particles;i_rank<SID.n_proc;i_rank++){
            n_groups_rank[i_rank]   =0;
            n_subgroups_rank[i_rank]=0;
            n_particles_rank[i_rank]=0;
            while(i_group<n_groups && n_particles_rank[i_rank]<(size_t)((float)(n_particles_left)/(float)(SID.n_proc-i_rank))){
              n_groups_rank[i_rank]++;
              n_subgroups_rank[i_rank]+=n_subgroups_group[i_group];
              n_particles_rank[i_rank]+=(size_t)group_length[i_group++];
            }
            n_particles_left-=(big_int)n_particles_rank[i_rank];
          }
        }
        else
          fseeko(fp,3*n_groups*sizeof(int),SEEK_CUR);

        // Report demcomposition results
        if(SID.n_proc>1){
          SID_log("Results of domain decomposition:",SID_LOG_OPEN);
          for(i_rank=0;i_rank<SID.n_proc;i_rank++) 
            SID_log("(i_rank,n_groups,n_subgroups,n_particles)=(%3d,%7d,%7d,%9lld)",SID_LOG_COMMENT,
                     i_rank,n_groups_rank[i_rank],n_subgroups_rank[i_rank],n_particles_rank[i_rank]);
          SID_log("Done.",SID_LOG_CLOSE);
        }

        // Send results to slave ranks
        SID_Bcast(&n_groups,                   sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(&n_subgroups,                sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(&n_particles,                sizeof(big_int),MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(group_length,     n_groups  *sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(group_offsets,    n_groups  *sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(n_subgroups_group,n_groups  *sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(n_groups_rank,    SID.n_proc*sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(n_subgroups_rank, SID.n_proc*sizeof(int),    MASTER_RANK,SID.COMM_WORLD);
        SID_Bcast(n_particles_rank, SID.n_proc*sizeof(size_t), MASTER_RANK,SID.COMM_WORLD);

        // Compute group offsets for this rank
        for(i_rank=0,rank_offset=0;i_rank<SID.My_rank;i_rank++) 
          rank_offset+=n_groups_rank[i_rank];

        // Store results
        ADaPS_store(&(plist->data),(void *)(&(group_length[rank_offset])),     "n_particles_group_%s",    ADaPS_COPY_SUBARRAY_INT,n_groups_rank[i_rank],catalog_name);
        ADaPS_store(&(plist->data),(void *)(&(group_offsets[rank_offset])),    "particle_offset_group_%s",ADaPS_COPY_SUBARRAY_INT,n_groups_rank[i_rank],catalog_name);
        ADaPS_store(&(plist->data),(void *)(&(n_subgroups_group[rank_offset])),"n_subgroups_group_%s",    ADaPS_COPY_SUBARRAY_INT,n_groups_rank[i_rank],catalog_name);
        SID_free((void **)(&group_length));     group_length     =(int *)ADaPS_fetch(plist->data,"n_particles_group_%s",    catalog_name);
        SID_free((void **)(&group_offsets));    group_offsets    =(int *)ADaPS_fetch(plist->data,"particle_offset_group_%s",catalog_name);
        SID_free((void **)(&n_subgroups_group));n_subgroups_group=(int *)ADaPS_fetch(plist->data,"n_subgroups_group_%s",    catalog_name);

        // Transform global offsets to local offsets
        for(i_group=0,offset_0=group_offsets[0];i_group<n_groups_rank[i_rank];i_group++)
          group_offsets[i_group]-=offset_0;

        SID_log("Done. (%d groups)",SID_LOG_CLOSE,n_groups);
      }
      else
        SID_log("NO GROUPS TO READ!",SID_LOG_CLOSE);
      // Store results
      ADaPS_store(&(plist->data),(void *)(&(n_groups_rank[SID.My_rank])),   "n_groups_%s",       ADaPS_SCALAR_INT,   catalog_name);
      ADaPS_store(&(plist->data),(void *)(&(n_groups)),                     "n_groups_all_%s",   ADaPS_SCALAR_INT,   catalog_name);
      ADaPS_store(&(plist->data),(void *)(&(n_subgroups_rank[SID.My_rank])),"n_subgroups_%s",    ADaPS_SCALAR_INT,   catalog_name);
      ADaPS_store(&(plist->data),(void *)(&(n_subgroups)),                  "n_subgroups_all_%s",ADaPS_SCALAR_INT,   catalog_name);
      ADaPS_store(&(plist->data),(void *)(&(n_particles_rank[SID.My_rank])),"n_particles_%s",    ADaPS_SCALAR_SIZE_T,catalog_name);
      ADaPS_store(&(plist->data),(void *)(&(n_particles)),                  "n_particles_all_%s",ADaPS_SCALAR_SIZE_T,catalog_name);
      fclose(fp);
    }
    else
      SID_trap_error("Could not open {%s}!\n",ERROR_LOGIC,filename_groups);
  }

  // Sub-group file
  if(flag_read_subgroups){
    SID_log("Reading subgroup catalogs {%s}...",SID_LOG_OPEN,filename_subgroups);
    if((fp=fopen(filename_subgroups,"r"))!=NULL){
      if(SID.I_am_Master)
        fread(&n_subgroups_in,sizeof(int),1,fp);
      else
        fseeko(fp,sizeof(int),SEEK_CUR);
      SID_Bcast(&n_subgroups_in,sizeof(int),MASTER_RANK,SID.COMM_WORLD);
      if(n_subgroups_in!=n_subgroups)
        SID_trap_error("n_subgroups_in(%d)!=n_subgroups(%d) in subgroup file {%s}!",ERROR_LOGIC,n_subgroups_in,n_subgroups,filename_subgroups);
      // Check that >0 groups are being loaded
      if(n_subgroups_in>0){
        // Allocate arrays
        subgroup_length=(int *)SID_malloc(sizeof(int)*n_subgroups_rank[SID.My_rank]);
        subgroup_offset=(int *)SID_malloc(sizeof(int)*n_subgroups_rank[SID.My_rank]);
        // Read subgroup lengths
        for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(i_rank==SID.My_rank)
            fread(subgroup_length,sizeof(int),n_subgroups_rank[i_rank],fp);
          else
            fseeko(fp,n_subgroups_rank[i_rank]*sizeof(int),SEEK_CUR);
        }
        // Read subgroup offsets
        for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(i_rank==SID.My_rank){
            fread(subgroup_offset,sizeof(int),n_subgroups_rank[i_rank],fp);
            // Transform global offsets to local offsets
            for(i_group=0,offset_0=subgroup_offset[0];i_group<n_subgroups_rank[i_rank];i_group++)
              subgroup_offset[i_group]-=offset_0;
          }
          else
            fseeko(fp,n_subgroups_rank[i_rank]*sizeof(int),SEEK_CUR);
        }
        // Count substructure particles
        for(i_group=0,n_particles=0;i_group<n_subgroups_rank[SID.My_rank];i_group++)
          n_particles+=subgroup_length[i_group];
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE,&n_particles,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
#endif
        // Store results
        ADaPS_store(&(plist->data),(void *)(subgroup_length),"n_particles_subgroup_%s",    ADaPS_DEFAULT,catalog_name);
        ADaPS_store(&(plist->data),(void *)(subgroup_offset),"particle_offset_subgroup_%s",ADaPS_DEFAULT,catalog_name);
        SID_log("Done. (%ld subgroups)",SID_LOG_CLOSE,n_subgroups_in);
      }
      else
        SID_log("NO SUBGROUPS TO READ!",SID_LOG_CLOSE);
      fclose(fp);
    }
    else
      SID_trap_error("Could not open {%s}!",ERROR_LOGIC,filename_subgroups);
  }
  
  // Read particle ids
  if(flag_read_ids){
    if((fp=fopen(filename_ids,"r"))!=NULL){
      SID_log("Reading group particle ids file {%s}...",SID_LOG_OPEN,filename_ids);
      fread(&id_byte_size, sizeof(int),1,fp);
      if(id_byte_size==sizeof(int)){
        fread(&n_ids_i,sizeof(int),1,fp);
        n_ids=(size_t)n_ids_i;
      }
      else
        fread(&n_ids,sizeof(size_t),1,fp);
      input_id=(size_t *)SID_malloc(sizeof(size_t)*n_particles_rank[SID.My_rank]);
      if(n_ids>0){
        // Read group particle lists
        for(i_rank=0;i_rank<SID.n_proc;i_rank++){
          if(i_rank==SID.My_rank){
            if(id_byte_size==sizeof(size_t))
              fread(input_id,id_byte_size,n_particles_rank[i_rank],fp);                        
            else{
              for(i_particle=0;i_particle<n_particles_rank[i_rank];i_particle++){
                fread(&id_in,id_byte_size,1,fp);
                input_id[i_particle]=(size_t)id_in;
              }              
            }
          }
          else
            fseeko(fp,n_particles_rank[i_rank]*id_byte_size,SEEK_CUR);
          SID_Barrier(SID.COMM_WORLD);
        }
        SID_log("Done. (%lld particles)",SID_LOG_CLOSE,n_ids);
      }
      else
        SID_log("NO IDS TO READ!",SID_LOG_CLOSE);
      // Store results
      ADaPS_store(&(plist->data),(void *)(input_id),"particle_ids_%s",ADaPS_DEFAULT,catalog_name);
      fclose(fp);
    }
    else
      SID_trap_error("Could not open {%s}!",ERROR_LOGIC,filename_ids);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);
  
  SID_profile_stop(SID_PROFILE_DEFAULT);
}
