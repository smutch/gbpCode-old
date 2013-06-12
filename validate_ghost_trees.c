#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpTrees.h>

int main(int argc,char *argv[]){
  SID_init(&argc,&argv,NULL);

  char filename_root[MAX_FILENAME_LENGTH];
  int  snap_start,snap_stop,snap_step;

  strcpy(filename_root,argv[1]);
  snap_start=atoi(argv[2]);
  snap_stop =atoi(argv[3]);
  snap_step =atoi(argv[4]);

  // Build ghost-halo version of horizontal tree files
  SID_log("Validating ghost trees...",SID_LOG_OPEN|SID_LOG_TIMER);

  // Loop over all the snapshots
  int snap_i;
  int i_snap;
  int n_groups_bad_total   =0;
  int n_subgroups_bad_total=0;
  for(snap_i=snap_stop,i_snap=0;snap_i>=snap_start;snap_i-=snap_step,i_snap++){
     FILE *fp_in;
     char  filename[MAX_FILENAME_LENGTH];
     SID_log("Processing snapshot No. %03d...",SID_LOG_OPEN,snap_i);

     // Open ghost group catalog file
     int i_file_group_ghosts;
     int n_file_group_ghosts;
     int n_group_ghosts;
     sprintf(filename,"%s/horizontal/ghost_catalogs/ghosts_%03d.catalog_subgroups_properties",filename_root,snap_i);
     if((fp_in=fopen(filename,"r"))==NULL)
        SID_trap_error("Could not open file {%s}",ERROR_IO_OPEN,filename);
     fread(&i_file_group_ghosts,sizeof(int),1,fp_in);
     fread(&n_file_group_ghosts,sizeof(int),1,fp_in);
     fread(&n_group_ghosts,     sizeof(int),1,fp_in);
     fclose(fp_in);

     // Open ghost group catalog file
     int i_file_subgroup_ghosts;
     int n_file_subgroup_ghosts;
     int n_subgroup_ghosts;
     sprintf(filename,"%s/horizontal/ghost_catalogs/ghosts_%03d.catalog_groups_properties",filename_root,snap_i);
     if((fp_in=fopen(filename,"r"))==NULL)
        SID_trap_error("Could not open file {%s}",ERROR_IO_OPEN,filename);
     fread(&i_file_subgroup_ghosts,sizeof(int),1,fp_in);
     fread(&n_file_subgroup_ghosts,sizeof(int),1,fp_in);
     fread(&n_subgroup_ghosts,     sizeof(int),1,fp_in);
     fclose(fp_in);

     // Open tree file
     sprintf(filename,"%s/horizontal/trees/horizontal_trees_ghosts_%03d.dat",filename_root,snap_i);
     if((fp_in=fopen(filename,"r"))==NULL)
        SID_trap_error("Could not open file {%s}",ERROR_IO_OPEN,filename);

     // Read tree header
     int n_step;
     int n_search;
     int n_groups;
     int n_subgroups;
     int n_groups_max;
     int n_subgroups_max; 
     int n_trees_subgroup;
     int n_trees_group;
     fread(&n_step,          sizeof(int),1,fp_in);
     fread(&n_search,        sizeof(int),1,fp_in);
     fread(&n_groups,        sizeof(int),1,fp_in);
     fread(&n_subgroups,     sizeof(int),1,fp_in);
     fread(&n_groups_max,    sizeof(int),1,fp_in);
     fread(&n_subgroups_max, sizeof(int),1,fp_in);
     fread(&n_trees_subgroup,sizeof(int),1,fp_in);
     fread(&n_trees_group,   sizeof(int),1,fp_in);

     // Check validity of indices
     int i_group;
     int i_subgroup;
     int n_group_descendant;
     int n_subgroup_descendant;
     int n_groups_bad   =0;
     int n_subgroups_bad=0;
     int n_null_ghost_groups=0;
     SID_log("No. of groups            = %d [%d ghost]",SID_LOG_COMMENT,n_groups,n_group_ghosts);
     SID_log("No. of subgroups         = %d [%d ghost]",SID_LOG_COMMENT,n_subgroups,n_subgroup_ghosts);
     for(i_group=0,i_subgroup=0;i_group<n_groups;i_group++){
  
        // Read groups
        int group_id;
        int group_type;
        int group_descendant_id;
        int group_tree_id;
        int group_file_offset;
        int group_index;
        int n_subgroups_group;
        fread(&group_id,           sizeof(int),1,fp_in);
        fread(&group_type,         sizeof(int),1,fp_in);
        fread(&group_descendant_id,sizeof(int),1,fp_in);
        fread(&group_tree_id,      sizeof(int),1,fp_in);
        //fread(&group_file_offset,  sizeof(int),1,fp_in);
        group_file_offset=1;
        fread(&group_index,        sizeof(int),1,fp_in);
        fread(&n_subgroups_group,  sizeof(int),1,fp_in);

        // Perform check on groups
        if(i_snap>0){
           if(group_index>=n_group_descendant){
             fprintf(stdout,"g: (%d,%d)->(%d,%d) %d\n",i_group,snap_i,group_index,snap_i+group_file_offset*snap_step,n_group_descendant);
             n_groups_bad++;
           }
        }

        if(check_mode_for_flag(group_type,TREE_CASE_GHOST_NULL))
           n_null_ghost_groups++;
       
        // Process subgroups 
        int j_subgroup;
        for(j_subgroup=0;j_subgroup<n_subgroups_group;i_subgroup++,j_subgroup++){
           if(i_subgroup>=n_subgroups)
              SID_trap_error("subgroup count exceeded at i_group=%d of %d",ERROR_LOGIC,i_group,n_groups);
  
           // Read subgroups
           int subgroup_id;
           int subgroup_type;
           int subgroup_descendant_id;
           int subgroup_tree_id;
           int subgroup_file_offset;
           int subgroup_index;
           fread(&subgroup_id,           sizeof(int),1,fp_in);
           fread(&subgroup_type,         sizeof(int),1,fp_in);
           fread(&subgroup_descendant_id,sizeof(int),1,fp_in);
           fread(&subgroup_tree_id,      sizeof(int),1,fp_in);
           //fread(&subgroup_file_offset,  sizeof(int),1,fp_in);
           subgroup_file_offset=1;
           fread(&subgroup_index,        sizeof(int),1,fp_in);

           // Perform check on subgroups
           if(i_snap>0){
              if(subgroup_index>=n_subgroup_descendant){
                fprintf(stdout,"s: (%d,%d)->(%d,%d) %d\n",i_subgroup,snap_i,subgroup_index,snap_i+subgroup_file_offset*snap_step,n_subgroup_descendant);
                n_subgroups_bad++;
              }
           }

           // Test message
           //if(check_mode_for_flag(group_type,TREE_CASE_GHOST) && snap_i==306)
           //   fprintf(stdout,"G: ( %d , %d , %d )->( %d , %d ) %d %d\n",
           //           i_group,i_subgroup,snap_i,subgroup_index,snap_i+subgroup_file_offset*snap_step,
           //           n_subgroup_descendant,subgroup_index>=n_subgroup_descendant);
        }

     }

     if(n_null_ghost_groups>0)
        SID_log("No. of null ghost groups = %d",SID_LOG_COMMENT,n_null_ghost_groups);

     if(n_groups_bad>0 || n_subgroups_bad>0){
        SID_log("No. of bad groups        = %d *****",SID_LOG_COMMENT,n_groups_bad);
        SID_log("No. of bad subgroups     = %d *****",SID_LOG_COMMENT,n_subgroups_bad);
     }
     n_groups_bad_total   +=n_groups_bad;
     n_subgroups_bad_total+=n_subgroups_bad;

     SID_log("Done.",SID_LOG_CLOSE);

     // Update descendant counts
     n_group_descendant   =n_groups;
     n_subgroup_descendant=n_subgroups;

     // Close file
     fclose(fp_in);
  }

  SID_log("Summary:",SID_LOG_OPEN);
  SID_log("No. of bad groups    = %d",SID_LOG_COMMENT,n_groups_bad_total);
  SID_log("No. of bad subgroups = %d",SID_LOG_COMMENT,n_subgroups_bad_total);
  SID_log("",SID_LOG_SILENT_CLOSE);

  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}     

