#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees.h>

void compute_trees_auxiliary(char *filename_root,
                             char *filename_snapshot_root,
                             char *filename_root_out,
                             int   i_read_start,
                             int   i_read_stop,
                             int   i_read_step,
                             int   n_search,
                             int   n_files_groups,
                             int   n_files_subgroups,
                             int  *flag_clean){
  int        i_read;
  int        i_snap;
  int        n_snap;
  int        k_match;
  int        n_write;
  int        i_write;
  int       *n_halos_snap;
  int       *halo_offset_snap;
  int       *n_halos_tree;
  int       *tree_MBP;
  int       *group_MBP;
  size_t    *group_MBP_index;
  size_t    *id_MBP_snap;
  size_t    *id_MBP_snap_index;
  int       *particle_halo;
  int       *tree_MBP_snap;
  size_t    *tree_MBP_snap_index;
  int        n_trees;
  int        n_halos;
  int        i_tree;
  int        j_tree;
  int        i_halo;
  int        j_halo;
  int        n_particles;
  int        n_list;
  int        n_ids;
  size_t     i_particle;
  size_t     j_particle;
  int        n_particles_snap;
  size_t     n_particles_mark_found;
  int        n_particles_snap_tree;
  int        offset_snap_tree;
  int        count_ids_halo;
  int        offset_ids_halo;
  size_t     min_id_halo;
  size_t     max_id_halo;
  int        min_index_halo;
  int        max_index_halo;
  SID_fp     fp_in;
  SID_fp     fp_out;
  SID_fp     fp_out_H;
  SID_fp     fp_out_A;
  SID_fp     fp_out_B;
  SID_fp     fp_out_C;
  SID_fp     fp_out_D;
  SID_fp     fp_out_E;
  SID_fp     fp_out_1;
  SID_fp     fp_out_2;
  SID_fp     fp_out_3;
  char       filename_in[256];
  char       filename_out[256];
  char       filename_snap[256];
  char       filename_temp_out_H[256];
  char       filename_temp_out_A[256];
  char       filename_temp_out_B[256];
  char       filename_temp_out_C[256];
  char       filename_temp_out_D[256];
  char       filename_temp_out_E[256];
  char       filename_temp_out_1[256];
  char       filename_temp_out_2[256];
  char       filename_temp_out_3[256];
  char       group_text_prefix[4];
  halo_MBP_info  halo_in;
  int            junk;
  size_t        *id_MBP;
  size_t        *id_MBP_index;
  GBPREAL          *x_MBP;
  GBPREAL          *y_MBP;
  GBPREAL          *z_MBP;
  GBPREAL          *vx_MBP;
  GBPREAL          *vy_MBP;
  GBPREAL          *vz_MBP;
  long long      id_temp;
  float          x_temp;
  float          y_temp;
  float          z_temp;
  float          vx_temp;
  float          vy_temp;
  float          vz_temp;

  SID_log("Constructing auxiliary tree files for snapshots #%d->#%d (step=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_start,i_read_stop,i_read_step);

  // Count how many snapshots we're using
  for(i_read=i_read_stop,n_snap=0;i_read>=i_read_start;i_read-=i_read_step) n_snap++;
  n_halos_snap    =(int *)SID_malloc(sizeof(int)*n_snap);
  halo_offset_snap=(int *)SID_malloc(sizeof(int)*n_snap);

  // Process subgroup trees and then group trees
  for(k_match=0;k_match<2;k_match++){
    switch(k_match){
      // Process subgroups
    case 0:
      n_write=n_files_subgroups;
      sprintf(group_text_prefix,"sub");
      break;
      // Process groups
    case 1:
      n_write=n_files_groups;
      sprintf(group_text_prefix,"");
      break;      
    }
    // Process each tree file in turn
    //   MODIFY FOR MPI: GIVE EACH CORE A FILE.
    SID_log("Writting %sgroup auxiliary files...",SID_LOG_OPEN|SID_LOG_TIMER,group_text_prefix);
    for(i_write=0;i_write<n_write;i_write++){

      // Open tree file and read header info
      if(n_write==1)
        sprintf(filename_in,"%s.%sgroup_trees_MBP",filename_root_out,group_text_prefix);
      else
        sprintf(filename_in,"%s.%sgroup_trees_MBP.%d",filename_root_out,group_text_prefix,i_write);
      SID_fopen(filename_in,"r",&fp_in);
      SID_fread(&n_trees,sizeof(int),1,&fp_in);
      SID_fread(&n_halos,sizeof(int),1,&fp_in);
      SID_log("Processing tree file {%s}...(n_trees=%d and n_halos=%d)...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in,n_trees,n_halos);

      // Allocate particle arrays and then finish reading header
      n_halos_tree      =(int    *)SID_malloc(sizeof(int)   *n_trees);
      id_MBP            =(size_t *)SID_malloc(sizeof(size_t)*n_halos);
      id_MBP_snap       =(size_t *)SID_malloc(sizeof(size_t)*n_halos);
      x_MBP             =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *n_halos);
      y_MBP             =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *n_halos);
      z_MBP             =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *n_halos);
      vx_MBP            =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *n_halos);
      vy_MBP            =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *n_halos);
      vz_MBP            =(GBPREAL   *)SID_malloc(sizeof(GBPREAL)  *n_halos);
      particle_halo     =(int    *)SID_malloc(sizeof(int)   *n_halos);
      tree_MBP          =(int    *)SID_malloc(sizeof(int)   *n_halos);
      group_MBP         =(int    *)SID_malloc(sizeof(int)   *n_halos);
      tree_MBP_snap     =(int    *)SID_malloc(sizeof(int)   *n_halos);

      SID_fread(n_halos_tree,sizeof(int),n_trees,&fp_in);

      // Read the halo list to count the number in each snapshot
      SID_log("Counting particles in each snapshot...",SID_LOG_OPEN|SID_LOG_TIMER);
      for(i_snap=0;i_snap<n_snap;i_snap++)
        n_halos_snap[i_snap]=0;
      for(j_tree=0,i_tree=0;j_tree<n_trees;j_tree++){
        for(i_halo=0;i_halo<n_halos_tree[j_tree];i_halo++,i_tree++){
          SID_fread(&halo_in,sizeof(halo_MBP_info),1,&fp_in);
          i_snap=((halo_in.snap_num)-i_read_start)/i_read_step;
          // Check that i_snap range is ok
          if(i_snap<0 || i_snap>=n_snap)  SID_trap_error("i_snap=%d! (%d)",ERROR_LOGIC,i_snap,n_snap);
          n_halos_snap[i_snap]++;
        }
      }
      SID_log("Done.",SID_LOG_CLOSE);
    
      // Set offsets
      halo_offset_snap[0]=0;
      for(i_snap=1;i_snap<n_snap;i_snap++)
        halo_offset_snap[i_snap]=halo_offset_snap[i_snap-1]+n_halos_snap[i_snap-1];

      // Reset file pointer to the start of the halo list
      SID_frewind(&fp_in);
      SID_fskip(sizeof(int),n_trees+2,&fp_in);

      // Now that arrays have been allocated and offsets determined, read the most-bount particle (MBP) ids, positions and velocities
      //   Place them in blocks ordered by snapshot
      SID_log("Reading particle IDs from the tree file...",SID_LOG_OPEN|SID_LOG_TIMER);
      for(i_snap=0;i_snap<n_snap;i_snap++)
        n_halos_snap[i_snap]=0;
      for(j_tree=0;j_tree<n_trees;j_tree++){
        for(i_halo=0;i_halo<n_halos_tree[j_tree];i_halo++){
          SID_fread(&halo_in,sizeof(halo_MBP_info),1,&fp_in);
          i_snap=((halo_in.snap_num)-i_read_start)/i_read_step;
          j_halo=halo_offset_snap[i_snap]+n_halos_snap[i_snap];
          // Sanity check: are i_snap and j_halo in the right range?
          if(i_snap<0 || i_snap>=n_snap)  SID_trap_error("i_snap=%d! (should be 0->%d)",ERROR_LOGIC,i_snap,n_snap-1);
          if(j_halo<0 || j_halo>=n_halos) SID_trap_error("j_halo=%d! (shoule be 0->%d)",ERROR_LOGIC,j_halo,n_halos-1);
          tree_MBP[j_halo] =j_tree;
          group_MBP[j_halo]=(size_t)(halo_in.group_halo_first);
          id_MBP[j_halo]   =(size_t)(halo_in.most_bound_id);
          x_MBP[j_halo]    =(size_t)(halo_in.pos[0]);
          y_MBP[j_halo]    =(size_t)(halo_in.pos[1]);
          z_MBP[j_halo]    =(size_t)(halo_in.pos[2]);
          vx_MBP[j_halo]   =(size_t)(halo_in.vel[0]);
          vy_MBP[j_halo]   =(size_t)(halo_in.vel[1]);
          vz_MBP[j_halo]   =(size_t)(halo_in.vel[2]);
          n_halos_snap[i_snap]++;
        }
      }
      SID_fclose(&fp_in);
      SID_log("Done.",SID_LOG_CLOSE);

      // Read MBP positions and velocities by processing each snapshot in turn; write to temporary files...
      SID_barrier(); // Make sure all ranks have read their tree file before all ranks begin reading all snapshots
      sprintf(filename_temp_out_H,"%s.%sgroup_hdH_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_A,"%s.%sgroup_hdA_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_B,"%s.%sgroup_hdB_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_C,"%s.%sgroup_hdC_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_D,"%s.%sgroup_hdD_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_E,"%s.%sgroup_hdE_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_1,"%s.%sgroup_ids_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_2,"%s.%sgroup_pos_temp.%d",filename_root_out,group_text_prefix,i_write);
      sprintf(filename_temp_out_3,"%s.%sgroup_vel_temp.%d",filename_root_out,group_text_prefix,i_write);
      SID_fopen(filename_temp_out_H,"w",&fp_out_H);
      SID_fopen(filename_temp_out_A,"w",&fp_out_A);
      SID_fopen(filename_temp_out_B,"w",&fp_out_B);
      SID_fopen(filename_temp_out_C,"w",&fp_out_C);
      SID_fopen(filename_temp_out_D,"w",&fp_out_D);
      SID_fopen(filename_temp_out_E,"w",&fp_out_E);
      SID_fopen(filename_temp_out_1,"w",&fp_out_1);
      SID_fopen(filename_temp_out_2,"w",&fp_out_2);
      SID_fopen(filename_temp_out_3,"w",&fp_out_3);
      
      // Write temporary id, position and velocity files
      for(i_snap=0,n_particles=0,n_ids=0,offset_snap_tree=0;i_snap<n_snap;i_snap++){
        SID_log("Processing snapshot #%d for auxilliary %sgroup tree file #%d...",SID_LOG_OPEN|SID_LOG_TIMER,i_read_stop-(n_snap-i_snap-1)*i_read_step,group_text_prefix,i_write+1);
        n_particles+=n_halos_snap[i_snap];

        // Write a junk array to the header (not sure why this is here!)
        junk=0;
        SID_fwrite(&junk,sizeof(int),1,&fp_out_A);
        SID_fwrite(&junk,sizeof(int),1,&fp_out_A);

        if(n_particles>0){
          // Create a sorted particle list
          SID_log("Initializing particle list...",SID_LOG_OPEN);
          merge_sort(id_MBP,(size_t)n_particles,&id_MBP_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          for(i_particle=0,j_particle=0;i_particle<n_particles;i_particle++){
            id_MBP_snap[j_particle]  =id_MBP[id_MBP_index[i_particle]];
            tree_MBP_snap[j_particle]=tree_MBP[id_MBP_index[i_particle]];
            j_particle++;
            if(j_particle>1){
              if(id_MBP_snap[j_particle-2]==id_MBP_snap[j_particle-1])
                j_particle--;
            }
          }
          n_particles_snap  =j_particle;
          n_ids            +=n_particles_snap;
          SID_log("%d particles...Done.",SID_LOG_CLOSE,n_particles_snap);

          // Write particle properties and count/offset arrays to temporary files -- snap_tree arrays first ...
          merge_sort(tree_MBP_snap,(size_t)n_particles_snap,&tree_MBP_snap_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          for(i_tree=0,i_particle=0;i_tree<n_trees;i_tree++){
            // Count the number of particles in the i_snap'th snapshot for the i_tree'th tree
            n_particles_snap_tree=0;
            if(n_particles_snap>0 && i_particle<n_particles_snap){
              while(tree_MBP_snap[tree_MBP_snap_index[i_particle]]<i_tree && i_particle<n_particles_snap-1) i_particle++;
              while(tree_MBP_snap[tree_MBP_snap_index[i_particle]]==i_tree && i_particle<n_particles_snap-1){
                n_particles_snap_tree++;
                i_particle++;
              }
              if(tree_MBP_snap[tree_MBP_snap_index[i_particle]]==i_tree){
                n_particles_snap_tree++;
                i_particle++;
              }
            }
            SID_fwrite(&n_particles_snap_tree,sizeof(int),1,&fp_out_B);
            SID_fwrite(&offset_snap_tree,     sizeof(int),1,&fp_out_C);
            offset_snap_tree+=n_particles_snap_tree;
          }

          // Sort particle IDs
          SID_log("Sort particle IDs...",SID_LOG_OPEN);
          merge_sort(id_MBP_snap,(size_t)n_particles_snap,&id_MBP_snap_index,SID_SIZE_T,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
          SID_log("Done.",SID_LOG_CLOSE);

          // Write arrays to temp files
          SID_log("Writing ids, positions and velocities to temporary files...",SID_LOG_OPEN|SID_LOG_TIMER);
          for(i_particle=0;i_particle<n_particles_snap;i_particle++){
            j_particle=id_MBP_snap_index[i_particle];
            id_temp=(long long)(id_MBP[j_particle]);
            x_temp =(float)(x_MBP[j_particle]);
            y_temp =(float)(y_MBP[j_particle]);
            z_temp =(float)(z_MBP[j_particle]);
            vx_temp=(float)(vx_MBP[j_particle]);
            vy_temp=(float)(vy_MBP[j_particle]);
            vz_temp=(float)(vz_MBP[j_particle]);
            SID_fwrite(&id_temp,sizeof(long long), 1,&fp_out_1);          
            SID_fwrite(&x_temp, sizeof(float),     1,&fp_out_2);          
            SID_fwrite(&y_temp, sizeof(float),     1,&fp_out_2);          
            SID_fwrite(&z_temp, sizeof(float),     1,&fp_out_2);
            SID_fwrite(&vx_temp,sizeof(float),     1,&fp_out_3);          
            SID_fwrite(&vy_temp,sizeof(float),     1,&fp_out_3);          
            SID_fwrite(&vz_temp,sizeof(float),     1,&fp_out_3);          
          }
          SID_log("Done.",SID_LOG_CLOSE);

          // Clean-up
          SID_free((void **)&id_MBP_snap_index);
          SID_free((void **)&tree_MBP_snap_index);
          SID_free((void **)&id_MBP_index);
        }
        // If there are no particles, write zeros to header arrays
        else{
          SID_log("NO HALOS...",SID_LOG_CONTINUE);
          junk=0;
          for(i_tree=0,i_particle=0;i_tree<n_trees;i_tree++){
            SID_fwrite(&junk,sizeof(int),1,&fp_out_B);
            SID_fwrite(&junk,sizeof(int),1,&fp_out_C);
          }    
        }
        
        SID_log("Done.",SID_LOG_CLOSE);
      } // i_snap

      // Write MBP counts and offsets for halo arrays ...
      SID_log("Writing halo arrays...",SID_LOG_OPEN|SID_LOG_TIMER);      
      merge_sort(group_MBP,(size_t)n_halos,&group_MBP_index,SID_INT,SORT_COMPUTE_INDEX,SORT_COMPUTE_NOT_INPLACE);
      for(i_halo=0,j_halo=0;i_halo<n_halos;i_halo++){

        // Find min and max MBP ids for this halo
        min_index_halo=n_halos-1;
        max_index_halo=0;
        while(i_halo==group_MBP[group_MBP_index[j_halo]] && j_halo<n_halos-1){
          min_index_halo=MIN(group_MBP_index[j_halo],min_index_halo);
          max_index_halo=MAX(group_MBP_index[j_halo],max_index_halo);
          j_halo++;
        }
        if(i_halo==group_MBP[group_MBP_index[j_halo]] && j_halo<n_halos-1){
          min_index_halo=MIN(group_MBP_index[j_halo],min_index_halo);
          max_index_halo=MAX(group_MBP_index[j_halo],max_index_halo);
          j_halo++;
        }

        // Convert to offsets and counts
        count_ids_halo=max_index_halo-min_index_halo+1;
        if(count_ids_halo>0)
          offset_ids_halo=min_index_halo;
        else{
          count_ids_halo = 0; // Without this, count gets set to -n_halos-2!
          offset_ids_halo=-1;
        }

        // Write results to temporary files
        SID_fwrite(&count_ids_halo,  sizeof(int),1,&fp_out_D);
        SID_fwrite(&offset_ids_halo, sizeof(int),1,&fp_out_E);
      }
      SID_free((void **)&group_MBP_index);
      SID_log("Done.",SID_LOG_CLOSE);

      // Top of header
      SID_fwrite(&n_halos,sizeof(int),1,&fp_out_H);          
      SID_fwrite(&n_ids,  sizeof(int),1,&fp_out_H);          
      SID_fwrite(&n_trees,sizeof(int),1,&fp_out_H);          
      SID_fwrite(&n_snap, sizeof(int),1,&fp_out_H);

      // Close temporary files
      SID_fclose(&fp_out_H);
      SID_fclose(&fp_out_A);
      SID_fclose(&fp_out_B);
      SID_fclose(&fp_out_C);
      SID_fclose(&fp_out_D);
      SID_fclose(&fp_out_E);
      SID_fclose(&fp_out_1);
      SID_fclose(&fp_out_2);
      SID_fclose(&fp_out_3);
      
      // Concatinate temporary files and write the final result    
      if(n_write==1)
        sprintf(filename_out,"%s.%sgroup_auxiliary",filename_root_out,group_text_prefix);
      else
        sprintf(filename_out,"%s.%sgroup_auxiliary.%d",filename_root_out,group_text_prefix,i_write);
      SID_cat_files(filename_out,
                    SID_CAT_CLEAN,9,
                    filename_temp_out_H,
                    filename_temp_out_A,
                    filename_temp_out_B,
                    filename_temp_out_C,
                    filename_temp_out_D,
                    filename_temp_out_E,
                    filename_temp_out_1,
                    filename_temp_out_2,
                    filename_temp_out_3);

      // Clean-up
      SID_free((void **)&n_halos_tree);
      SID_free((void **)&id_MBP);
      SID_free((void **)&x_MBP);
      SID_free((void **)&y_MBP);
      SID_free((void **)&z_MBP);
      SID_free((void **)&vx_MBP);
      SID_free((void **)&vy_MBP);
      SID_free((void **)&vz_MBP);
      SID_free((void **)&particle_halo);
      SID_free((void **)&tree_MBP);
      SID_free((void **)&tree_MBP_snap);

      SID_log("Done.",SID_LOG_CLOSE);
    } // i_write
    SID_log("Done.",SID_LOG_CLOSE);

  } // k_match
  
  // Clean-up
  SID_free((void **)&n_halos_snap);
  SID_free((void **)&halo_offset_snap);
  SID_log("Done.",SID_LOG_CLOSE);

}
