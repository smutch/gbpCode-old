#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>
#include <gbpHalos.h>

int main(int argc, char *argv[]){
  plist_info  plist;
  char        filename_root[256];
  char        filename_log[256];
  char        filename_number[256];
  char        filename_in_halos[256];
  char        filename_out_groups[256];
  char        filename_out_groups_A[256];
  char        filename_out_groups_B[256];
  char        filename_out_groups_C[256];
  char        filename_out_subgroups[256];
  char        filename_out_subgroups_A[256];
  char        filename_out_subgroups_B[256];
  char        filename_out_hierarchy[256];
  char        filename_out_hierarchy_A[256];
  char        filename_out_hierarchy_B[256];
  char        filename_out_particles[256];
  char        i_match_txt[5];
  int         n_groups_AHF;
  int         n_groups;
  int         n_subgroups;
  int         n_subgroups_matched;
  int         n_subgroups_group;
  size_t      n_particles;
  size_t      n_particles_in_groups;
  size_t      n_particles_in_subgroups;
  size_t      n_particles_AHF_not_used;
  int         n_particles_temp;
  int        *n_p_1=NULL;
  int         flag_continue;
  int         flag_long_ids;
  int         i_match;
  int         match_id_next;
  int        *match_id=NULL;
  int        *match_id_initial=NULL;
  FILE       *fp=NULL;
  FILE       *fp_in_halos=NULL;
  FILE       *fp_out=NULL;
  int         n_match;
  int        *id_2=NULL;
  id_int     *particle_ids_AHF=NULL;
  size_t     *particle_ids_AHF_index=NULL;
  id_int      id_largest;
  int         id_byte_size;
  size_t     *group_particles=NULL;
  int         group_id;
  int         subgroup_id;
  int         i_group;
  int         j_group;
  int         k_group;
  size_t     n_particles_AHF;
  int        *subgroup_size=NULL;
  int        *hierarchy_level=NULL;
  int        *hierarchy_match=NULL;
  int         subgroup_size_max;
  int        *subgroup_size_list=NULL;
  int        *subgroup_index_list=NULL;
  size_t     *subgroup_size_list_index=NULL;
  int        *group_offsets=NULL;
  size_t      group_index;
  int        *group_size=NULL;
  int        *group_size_AHF=NULL;
  int        *group_offsets_AHF=NULL;
  int         max_subgroup_size;
  int         i_subgroup;
  int         j_subgroup;
  int         n_subgroups_group_max;
  size_t     *group_size_index=NULL;
  size_t     *match_id_index=NULL;
  size_t      subgroup_index;
  int         group_offset;
  int         subgroup_offset;
  int         group_count;
  size_t     *group_particles_index=NULL;
  size_t     *subgroup_particles=NULL;
  int        *particle_group=NULL;
  size_t     *particle_group_index=NULL;
  size_t      i_particle;
  size_t      j_particle;
  size_t      k_particle;
  int         i_file;
  int         i_file_start;
  int         i_file_stop;
  size_t     *match_index=NULL;
  int         flag_match_subgroups;
  FILE       *fp_log=NULL;
  FILE       *fp_in=NULL;
  FILE       *fp_out_particles=NULL;
  FILE       *fp_out_groups=NULL;
  FILE       *fp_out_groups_A=NULL;
  FILE       *fp_out_groups_B=NULL;
  FILE       *fp_out_groups_C=NULL;
  FILE       *fp_out_subgroups_A=NULL;
  FILE       *fp_out_subgroups_B=NULL;
  FILE       *fp_out_hierarchy_A=NULL;
  FILE       *fp_out_hierarchy_B=NULL;
  FILE       *fp_test=NULL;
  int         substructure_level;
  int         substructure_level_max;
  halo_properties_info *properties=NULL;
  void                 *particle_buffer=NULL;
  int                   flag_found;

  SID_init(&argc,&argv,NULL);

  strcpy(filename_root,argv[1]);
  i_file_start=atoi(argv[2]);
  i_file_stop =atoi(argv[3]);
  
  SID_log("Converting files #%d->#%d from AHF to subfind format...",
    SID_LOG_OPEN|SID_LOG_TIMER,
    i_file_start,i_file_stop);

  sprintf(filename_log,"%s_%dto%d.convert_AHF_log",filename_root,i_file_start,i_file_stop);

  // Loop over all files
  for(i_file=i_file_start;i_file<=i_file_stop;i_file++){
    SID_log("Processing file #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_file);

    // Read catalogs
    if(i_file<10)
      sprintf(filename_number,"00%1d",i_file);
    else if(i_file<100)
      sprintf(filename_number,"0%2d", i_file);
    else
      sprintf(filename_number,"%3d", i_file);

    // Read AHF group file
    init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);
    read_groups_AHF(filename_root,i_file,READ_GROUPS_ALL,&plist,filename_number);
    n_groups_AHF     =((int     *)ADaPS_fetch(plist.data,"n_groups_%s",   filename_number))[0];

    n_groups                =0;
    n_subgroups             =0;
    n_subgroups_matched     =0;
    n_subgroups_group       =0;
    n_particles             =0;
    n_particles_in_groups   =0;
    n_particles_in_subgroups=0;
    n_particles_AHF_not_used=0;
    n_subgroups_group_max   =0;

    if(n_groups_AHF>0){
    n_particles_AHF  =(size_t)((size_t *)ADaPS_fetch(plist.data,"n_particles_%s",filename_number))[0];
    group_size_AHF   =(int    *)ADaPS_fetch(plist.data,"n_particles_group_%s",     filename_number);
    group_offsets_AHF=(int    *)ADaPS_fetch(plist.data,"particle_offset_group_%s", filename_number);
    particle_ids_AHF =(size_t *)ADaPS_fetch(plist.data,"particle_ids_%s",          filename_number);

    // Find largest id so we know what size to write the ids with
    for(i_particle=0,id_largest=0;i_particle<n_particles_AHF;i_particle++)
      id_largest=MAX(id_largest,particle_ids_AHF[i_particle]);
    if(id_largest>INT_MAX){
      flag_long_ids=TRUE;
      id_byte_size =sizeof(size_t);
    }
    else{
      flag_long_ids=FALSE;
      id_byte_size =sizeof(int);
    }

    // Match AHF groups against themselves to find substructure
    match_halos(&plist,i_file,NULL,0,&plist,i_file,NULL,0,"substructure",MATCH_SUBSTRUCTURE);
    match_id_initial=(int  *)ADaPS_fetch(plist.data,"match_substructure");
    hierarchy_match =match_id_initial; // Fore readability

    // Assign sub-...-sub-structures to parent (ie. top-level) halos
    SID_log("Assigning substructures to groups...",SID_LOG_OPEN);
    group_size     =(int *)SID_malloc(sizeof(int)*n_groups_AHF);
    subgroup_size  =(int *)SID_malloc(sizeof(int)*n_groups_AHF);
    hierarchy_level=(int *)SID_malloc(sizeof(int)*n_groups_AHF);
    particle_group =(int *)SID_malloc(sizeof(int)*n_particles_AHF);
    for(i_group=0,i_particle=0;i_group<n_groups_AHF;i_group++){
      group_size[i_group]   =0;
      subgroup_size[i_group]=0;
      for(j_particle=0;j_particle<group_size_AHF[i_group];i_particle++,j_particle++)
        particle_group[i_particle]=i_group;
    }
    match_id=(int *)SID_malloc(sizeof(int)*n_groups_AHF);
    for(i_group=0,substructure_level_max=0;i_group<n_groups_AHF;i_group++){
      substructure_level=0;
      match_id_next     =match_id_initial[i_group];
      match_id[i_group] =match_id_next;
      while(match_id_next>=0){
        substructure_level++;
        match_id[i_group]=match_id_next; // Tie subgroups to their top-level group
        match_id_next    =match_id_initial[match_id_next];
      }
      if(match_id[i_group]<0) 
        match_id[i_group]=i_group; // Unmatched halos should be matched to themselves
      hierarchy_level[i_group]=substructure_level;
      substructure_level_max=MAX(substructure_level,substructure_level_max);
    }
    // needed? ADaPS_store(&(plist.data),(void *)(match_id),"match_substructure",ADaPS_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);

    // Make sure the deepest substructures are given particle ownership
    SID_log("Assigning particles to subgroups...",SID_LOG_OPEN);
    merge_sort((void *)particle_ids_AHF,(size_t)n_particles_AHF,&particle_ids_AHF_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    for(i_particle=0,n_particles_AHF_not_used=0;i_particle<n_particles_AHF;i_particle+=k_particle){
      // Count the number of times this particle id is used
      j_particle=i_particle;
      while(particle_ids_AHF[particle_ids_AHF_index[j_particle]]==particle_ids_AHF[particle_ids_AHF_index[i_particle]] && j_particle<(n_particles_AHF-2)) j_particle++;
      if(particle_ids_AHF[particle_ids_AHF_index[j_particle]]==particle_ids_AHF[particle_ids_AHF_index[i_particle]])                                      j_particle++;
      k_particle=j_particle-i_particle;
      // Find the deepest substructure using this particle id...
      i_group=particle_group[particle_ids_AHF_index[i_particle]];
      for(j_particle=1;j_particle<k_particle;j_particle++){
        j_group=particle_group[particle_ids_AHF_index[i_particle+j_particle]];
        if(group_size_AHF[j_group]<group_size_AHF[i_group])
          i_group=j_group;
      }
      // ... and set particle's group to a dummy value if this particle instance is not from the deepest group
      for(j_particle=0,flag_found=FALSE;j_particle<k_particle;j_particle++){
        if(particle_group[particle_ids_AHF_index[i_particle+j_particle]]!=i_group || flag_found){
          particle_group[particle_ids_AHF_index[i_particle+j_particle]]=-1;
          n_particles_AHF_not_used++;
        }
        else
          flag_found=TRUE;
      }
    }
    SID_free((void **)&particle_ids_AHF_index);
    SID_log("Done.",SID_LOG_CLOSE);

    // Generate subgroup_size array
    for(i_group=0;i_group<n_groups_AHF;i_group++)
      subgroup_size[i_group]=0;
    for(i_particle=0;i_particle<n_particles_AHF;i_particle++){
      i_group=particle_group[i_particle];
      if(i_group>=0)
        subgroup_size[i_group]++;
    }

    // Get rid of groups that are too small
    for(i_particle=0;i_particle<n_particles_AHF;i_particle++){
      i_group=particle_group[i_particle];
      if(i_group>=0){
        if(subgroup_size[i_group]<20){
          n_particles_AHF_not_used++;
          particle_group[i_particle]=-1;
        }
      }
    }

    // Regenerate subgroup_size array
    for(i_group=0;i_group<n_groups_AHF;i_group++)
      subgroup_size[i_group]=0;
    for(i_particle=0;i_particle<n_particles_AHF;i_particle++){
      i_group=particle_group[i_particle];
      if(i_group>=0)
        subgroup_size[i_group]++;
    }

    // Find the largest subgroup's size
    for(i_group=0,n_subgroups=0,subgroup_size_max=0;i_group<n_groups_AHF;i_group++)
      subgroup_size_max=MAX(subgroup_size[i_group],subgroup_size_max);

    // Generate group_size array
    for(i_group=0;i_group<n_groups_AHF;i_group++)
      group_size[match_id[i_group]]+=subgroup_size[i_group]; // update group size

    // Sort groups in order of size
    merge_sort((void *)group_size,(size_t)n_groups_AHF,&group_size_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
    merge_sort((void *)match_id,  (size_t)n_groups_AHF,&match_id_index,  SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

    // Count groups, subgroups, etc.
    SID_log("Counting groups & subgroups...",SID_LOG_OPEN);
    for(i_group=0,n_groups=0,n_subgroups=0;i_group<n_groups_AHF;i_group++){
      group_index =group_size_index[n_groups_AHF-i_group-1];
      
      // Find start of subgroup list for this group
      j_group=find_index_int(match_id,group_index,n_groups_AHF,match_id_index);
      while(group_index>match_id[match_id_index[j_group]] && j_group<(n_groups_AHF-2)) j_group++;
      if(group_index>match_id[match_id_index[j_group]])                                j_group++;

      // Count subgroups
      n_subgroups_group=0;
      while(match_id[match_id_index[j_group]]==group_index && j_group<(n_groups_AHF-2)){
        if(subgroup_size[match_id_index[j_group]]>0)
          n_subgroups_group++;
        j_group++;
      }
      if(match_id[match_id_index[j_group]]==group_index){
        if(subgroup_size[match_id_index[j_group]]>0)
          n_subgroups_group++;
        j_group++;        
      }
      n_subgroups+=n_subgroups_group;

      // Largest number of subgroups
      n_subgroups_group_max=MAX(n_subgroups_group_max,n_subgroups_group);    

      // Count groups
      if(n_subgroups_group>0)
        n_groups++;
    }
    SID_log("Done.",SID_LOG_CLOSE);
    }

    // Find largest subgroup and count the number of particles in groups
    for(i_group=0,max_subgroup_size=0,n_particles_in_groups=0;i_group<n_groups_AHF;i_group++){
      max_subgroup_size=MAX(max_subgroup_size,subgroup_size[i_group]);
      if(subgroup_size[i_group]>0)
        n_particles_in_groups+=(size_t)subgroup_size[i_group];
    }

    // Write some statistics
    SID_log("Substructure statistics:",SID_LOG_OPEN);
    SID_log("Number of groups                 =%d",          SID_LOG_COMMENT,n_groups);
    SID_log("Number of subgroups              =%d",          SID_LOG_COMMENT,n_subgroups);
    SID_log("Max number of subgroups per group=%d",          SID_LOG_COMMENT,n_subgroups_group_max);
    SID_log("Largest subgroup                 =%d particles",SID_LOG_COMMENT,subgroup_size_max);
    SID_log("Depth of substructure heirarchy  =%d levels",   SID_LOG_COMMENT,substructure_level_max);
    SID_log("Number of AHF particles used     =%lld",        SID_LOG_COMMENT,n_particles_in_groups);
    SID_log("Number of AHF particles NOT used =%lld",        SID_LOG_COMMENT,n_particles_AHF_not_used);
    SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);     

    // Open files
    SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,0);
    SID_log("Writing %d groups, %d subgroups and %lld particles to files...",SID_LOG_OPEN|SID_LOG_TIMER,n_groups,n_subgroups,n_particles_in_groups);
    sprintf(filename_out_groups,     "%s_%s.catalog_groups",     filename_root,filename_number);
    sprintf(filename_out_groups_A,   "%s_%s.catalog_groups_A",   filename_root,filename_number);
    sprintf(filename_out_groups_B,   "%s_%s.catalog_groups_B",   filename_root,filename_number);
    sprintf(filename_out_groups_C,   "%s_%s.catalog_groups_C",   filename_root,filename_number);
    sprintf(filename_out_subgroups,  "%s_%s.catalog_subgroups",  filename_root,filename_number);
    sprintf(filename_out_subgroups_A,"%s_%s.catalog_subgroups_A",filename_root,filename_number);
    sprintf(filename_out_subgroups_B,"%s_%s.catalog_subgroups_B",filename_root,filename_number);
    sprintf(filename_out_hierarchy,  "%s_%s.catalog_hierarchy",  filename_root,filename_number);
    sprintf(filename_out_hierarchy_A,"%s_%s.catalog_hierarchy_A",filename_root,filename_number);
    sprintf(filename_out_hierarchy_B,"%s_%s.catalog_hierarchy_B",filename_root,filename_number);
    sprintf(filename_out_particles,  "%s_%s.catalog_particles",  filename_root,filename_number);
    fp_out_groups_A   =fopen(filename_out_groups_A,   "w");
    fp_out_groups_B   =fopen(filename_out_groups_B,   "w");
    fp_out_groups_C   =fopen(filename_out_groups_C,   "w");
    fp_out_subgroups_A=fopen(filename_out_subgroups_A,"w");
    fp_out_subgroups_B=fopen(filename_out_subgroups_B,"w");
    fp_out_hierarchy_A=fopen(filename_out_hierarchy_A,"w");
    fp_out_hierarchy_B=fopen(filename_out_hierarchy_B,"w");
    fp_out_particles  =fopen(filename_out_particles,  "w");

    // Write headers
    fwrite(&n_groups,              sizeof(int),    1,fp_out_groups_A);
    fwrite(&n_subgroups,           sizeof(int),    1,fp_out_subgroups_A);
    fwrite(&n_subgroups,           sizeof(int),    1,fp_out_hierarchy_A);
    fwrite(&id_byte_size,          sizeof(int),    1,fp_out_particles);
    switch(flag_long_ids){
      case TRUE:
        fwrite(&n_particles_in_groups, sizeof(size_t),1,fp_out_particles);
        break;
      default:
        n_particles_temp=(int)n_particles_in_groups;
        fwrite(&n_particles_temp,sizeof(int),1,fp_out_particles);
        break;
    }

    // Write files; group and subgroup files in parts (to be concatinated together later)
    subgroup_size_list =(int  *)SID_malloc(sizeof(int) *n_subgroups_group_max);
    subgroup_index_list=(int  *)SID_malloc(sizeof(int) *n_subgroups_group_max);
    particle_buffer    =(void *)SID_malloc(id_byte_size*subgroup_size_max);
    subgroup_offset    =0;
    group_offset       =0;

    for(i_group=n_groups_AHF-1;i_group>=n_groups_AHF-n_groups;i_group--){
      group_index=group_size_index[i_group];

      // Find start of subgroup list for this group
      i_subgroup=find_index_int(match_id,group_index,n_groups_AHF,match_id_index);
      while(group_index>match_id[match_id_index[i_subgroup]] && i_subgroup<(n_groups_AHF-2)) i_subgroup++;
      if(group_index>match_id[match_id_index[i_subgroup]])                                   i_subgroup++;

      // Create a list of subgroups for this group and sort it by size
      n_subgroups_group=0;
      subgroup_index   =match_id_index[i_subgroup];
      while(match_id[subgroup_index]==group_index && i_subgroup<(n_groups_AHF-2)){
        if(subgroup_size[subgroup_index]>0){
          subgroup_size_list[n_subgroups_group] =subgroup_size[subgroup_index];
          subgroup_index_list[n_subgroups_group]=(int)subgroup_index;
          n_subgroups_group++;
        }
        i_subgroup++;
        subgroup_index=match_id_index[i_subgroup];
      }
      if(match_id[subgroup_index]==group_index){
        if(subgroup_size[subgroup_index]>0){
          subgroup_size_list[n_subgroups_group] =subgroup_size[subgroup_index];
          subgroup_index_list[n_subgroups_group]=(int)subgroup_index;
          n_subgroups_group++;
        }
        i_subgroup++;
        subgroup_index=match_id_index[i_subgroup];
      }
      merge_sort((void *)subgroup_size_list,(size_t)n_subgroups_group,&subgroup_size_list_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);

      // Perform writes for subgroups and particle lists
      for(i_subgroup=0,i_particle=0;i_subgroup<n_subgroups_group;i_subgroup++){
        j_subgroup=subgroup_index_list[subgroup_size_list_index[n_subgroups_group-i_subgroup-1]];
        // ... subgroups ...
        fwrite(&(subgroup_size[j_subgroup]),   sizeof(int),1,fp_out_subgroups_A);
        fwrite(&(subgroup_offset),             sizeof(int),1,fp_out_subgroups_B);
        fwrite(&(hierarchy_match[j_subgroup]), sizeof(int),1,fp_out_hierarchy_A);
        fwrite(&(hierarchy_level[j_subgroup]), sizeof(int),1,fp_out_hierarchy_B);
        subgroup_offset+=subgroup_size[j_subgroup];
        // ... and particles
        for(j_particle=group_offsets_AHF[j_subgroup],k_particle=0,i_particle=0;k_particle<group_size_AHF[j_subgroup];j_particle++,k_particle++){
          if(particle_group[j_particle]==j_subgroup){
            switch(flag_long_ids){
              case TRUE:
                ((size_t *)particle_buffer)[i_particle++]=(size_t)(particle_ids_AHF[j_particle]);
                break;
              default:
                ((int *)particle_buffer)[i_particle++]=(int)(particle_ids_AHF[j_particle]);
                break;
            }
          }
        }
        if(i_particle==subgroup_size[j_subgroup])
          fwrite(particle_buffer,id_byte_size,i_particle,fp_out_particles);
        else
          SID_trap_error("Subgroup size mismatch!",ERROR_LOGIC);
      }

      SID_free((void **)&subgroup_size_list_index);

      // Perform writes for groups
      fwrite(&(group_size[group_index]), sizeof(int),1,fp_out_groups_A);
      fwrite(&group_offset,              sizeof(int),1,fp_out_groups_B);
      fwrite(&n_subgroups_group,         sizeof(int),1,fp_out_groups_C);
      group_offset+=group_size[group_index];
    }
    SID_free((void **)&subgroup_size_list);    
    SID_free((void **)&subgroup_index_list);    
    SID_free((void **)&particle_buffer);    

    fclose(fp_out_groups_A);
    fclose(fp_out_groups_B);
    fclose(fp_out_groups_C);
    fclose(fp_out_subgroups_A);
    fclose(fp_out_subgroups_B);
    fclose(fp_out_hierarchy_A);
    fclose(fp_out_hierarchy_B);
    fclose(fp_out_particles);

    // Concatinate group and subgroup temp files into final files
    SID_cat_files(filename_out_groups,SID_CAT_CLEAN,3,
                  filename_out_groups_A,
                  filename_out_groups_B,
                  filename_out_groups_C);

    SID_cat_files(filename_out_subgroups,SID_CAT_CLEAN,2,
                  filename_out_subgroups_A,
                  filename_out_subgroups_B);

    SID_cat_files(filename_out_hierarchy,SID_CAT_CLEAN,2,
                  filename_out_hierarchy_A,
                  filename_out_hierarchy_B);

    // Clean-up
    SID_free((void **)&subgroup_size);
    SID_free((void **)&hierarchy_level);
    SID_free((void **)&group_size);
    SID_free((void **)&group_size_index);
    SID_free((void **)&match_id_index);
    free_plist(&plist);
    SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
    SID_log("Done.",SID_LOG_CLOSE);

    // Write log file
    SID_log("Writing to log file...",SID_LOG_OPEN);
    // Write a header for the log file
    if(i_file==i_file_start){
      fp_log=fopen(filename_log,"w");
      fprintf(fp_log,"# (1):  filenumber\n");
      fprintf(fp_log,"# (2):  n_groups_AHF\n");
      fprintf(fp_log,"# (3):  n_particles_AHF\n");
      fprintf(fp_log,"# (4):  n_groups\n");
      fprintf(fp_log,"# (5):  n_subgroups\n");
      fprintf(fp_log,"# (6):  max number of subgroups per group\n");
      fprintf(fp_log,"# (7):  largest subgroup\n");
      fprintf(fp_log,"# (8):  depth of substructure heirarchy\n");      
      fprintf(fp_log,"# (9):  number of AHF particles used\n");      
      fprintf(fp_log,"# (10): number of AHF particles NOT used\n");      
    }
    else
      fp_log=fopen(filename_log,"a");
    fprintf(fp_log,"%4d %9d %12zd %9d %9d %9d %9d %9d %12zd %12zd\n",
      i_file,
      n_groups_AHF,
      n_particles_AHF,
      n_groups,
      n_subgroups,
      n_subgroups_group_max,
      subgroup_size_max,
      substructure_level_max,
      n_particles_in_groups,
      n_particles_AHF_not_used);
    fclose(fp_log);
    
    SID_log("Done.",SID_LOG_CLOSE);    

    SID_log("Done.",SID_LOG_CLOSE);
  }
  
  SID_log("Done.",SID_LOG_CLOSE);  
  SID_exit(ERROR_NONE);
}
