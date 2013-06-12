#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpTrees.h>

void interpolate_halo_local(halo_properties_info *halo_initial,
                            halo_properties_info *halo_final,
                            halo_properties_info *halo_return,
                            double                time_initial,
                            double                time_final,
                            double                time_return,int i_group,int mode,int offset,const char *txt);
void interpolate_halo_local(halo_properties_info *halo_initial,
                            halo_properties_info *halo_final,
                            halo_properties_info *halo_return,
                            double                time_initial,
                            double                time_final,
                            double                time_return,int i_group,int mode,int offset,const char *txt){

   if(halo_initial==NULL)
      SID_trap_error("Initial halo is undefined in interpolate_halo_local(). i_%sgroup=%d mode=%d offset=%d",ERROR_LOGIC,txt,i_group,mode,offset);
   if(halo_final==NULL)
      SID_trap_error("Final halo is undefined in interpolate_halo_local(). i_%sgroup=%d mode=%d offset=%d",ERROR_LOGIC,txt,i_group,mode,offset);

   double f_interpolate;
   f_interpolate=(time_return-time_initial)/(time_final-time_initial);

   // Pass the MBP forward
   halo_return->id_MBP=halo_initial->id_MBP;

   // Perform interpolation
   halo_return->M_vir          =halo_initial->M_vir+f_interpolate*(halo_final->M_vir-halo_initial->M_vir);
   halo_return->n_particles    =0;
   halo_return->position_COM[0]=halo_initial->position_COM[0]+f_interpolate*(halo_final->position_COM[0]-halo_initial->position_COM[0]);
   halo_return->position_COM[1]=halo_initial->position_COM[1]+f_interpolate*(halo_final->position_COM[1]-halo_initial->position_COM[1]);
   halo_return->position_COM[2]=halo_initial->position_COM[2]+f_interpolate*(halo_final->position_COM[2]-halo_initial->position_COM[2]);
   halo_return->velocity_COM[0]=halo_initial->velocity_COM[0]+f_interpolate*(halo_final->velocity_COM[0]-halo_initial->velocity_COM[0]);
   halo_return->velocity_COM[1]=halo_initial->velocity_COM[1]+f_interpolate*(halo_final->velocity_COM[1]-halo_initial->velocity_COM[1]);
   halo_return->velocity_COM[2]=halo_initial->velocity_COM[2]+f_interpolate*(halo_final->velocity_COM[2]-halo_initial->velocity_COM[2]);
   halo_return->position_MBP[0]=halo_initial->position_MBP[0]+f_interpolate*(halo_final->position_MBP[0]-halo_initial->position_MBP[0]);
   halo_return->position_MBP[1]=halo_initial->position_MBP[1]+f_interpolate*(halo_final->position_MBP[1]-halo_initial->position_MBP[1]);
   halo_return->position_MBP[2]=halo_initial->position_MBP[2]+f_interpolate*(halo_final->position_MBP[2]-halo_initial->position_MBP[2]);
   halo_return->velocity_MBP[0]=halo_initial->velocity_MBP[0]+f_interpolate*(halo_final->velocity_MBP[0]-halo_initial->velocity_MBP[0]);
   halo_return->velocity_MBP[1]=halo_initial->velocity_MBP[1]+f_interpolate*(halo_final->velocity_MBP[1]-halo_initial->velocity_MBP[1]);
   halo_return->velocity_MBP[2]=halo_initial->velocity_MBP[2]+f_interpolate*(halo_final->velocity_MBP[2]-halo_initial->velocity_MBP[2]);
   halo_return->spin[0]        =halo_initial->spin[0]        +f_interpolate*(halo_final->spin[0]        -halo_initial->spin[0]);
   halo_return->spin[1]        =halo_initial->spin[1]        +f_interpolate*(halo_final->spin[1]        -halo_initial->spin[1]);
   halo_return->spin[2]        =halo_initial->spin[2]        +f_interpolate*(halo_final->spin[2]        -halo_initial->spin[2]);
   halo_return->R_vir          =halo_initial->R_vir          +f_interpolate*(halo_final->R_vir          -halo_initial->R_vir);
   halo_return->R_halo         =halo_initial->R_halo         +f_interpolate*(halo_final->R_halo         -halo_initial->R_halo);
   halo_return->R_max          =halo_initial->R_max          +f_interpolate*(halo_final->R_max          -halo_initial->R_max);
   halo_return->V_max          =halo_initial->V_max          +f_interpolate*(halo_final->V_max          -halo_initial->V_max);
   halo_return->sigma_v        =halo_initial->sigma_v        +f_interpolate*(halo_final->sigma_v        -halo_initial->sigma_v);
   halo_return->q_triaxial     =halo_initial->q_triaxial     +f_interpolate*(halo_final->q_triaxial     -halo_initial->q_triaxial);
   halo_return->s_triaxial     =halo_initial->s_triaxial     +f_interpolate*(halo_final->s_triaxial     -halo_initial->s_triaxial);
   int i_vec;
   int j_vec;
   for(i_vec=0;i_vec<3;i_vec++){
      for(j_vec=0;j_vec<3;j_vec++){
         halo_return->shape_eigen_vectors[i_vec][j_vec]=
            halo_initial->shape_eigen_vectors[i_vec][j_vec]+
            f_interpolate*(halo_final->shape_eigen_vectors[i_vec][j_vec]-halo_initial->shape_eigen_vectors[i_vec][j_vec]);
      }
   }
}

void write_ghost_catalog(tree_horizontal_ghost_group_info      *groups,
                         halo_properties_info                ***group_properties,
                         int                                    n_group_ghosts,
                         int                                    n_groups,
                         tree_horizontal_ghost_subgroup_info   *subgroups,
                         halo_properties_info                ***subgroup_properties,
                         int                                    n_subgroup_ghosts,
                         int                                    n_subgroups,
                         char                                  *filename_output_dir,
                         char                                  *filename_cat,
                         int                                    i_write,
                         int                                    j_write,
                         int                                    l_write,
                         int                                    n_wrap,
                         double                                *a_list,
                         cosmo_info                           **cosmo){

   int  i_loop;
   halo_properties_info  *group_return;
   halo_properties_info  *subgroup_return;
   group_return   =(halo_properties_info *)SID_calloc(sizeof(halo_properties_info));
   subgroup_return=(halo_properties_info *)SID_calloc(sizeof(halo_properties_info));
   SID_log("Writing ghost catalogs for snapshot #%03d...",SID_LOG_OPEN,j_write);

   // Open files and write headers
   SID_fp fp_groups;
   SID_fp fp_subgroups;
   char   filename_output_dir_horizontal[MAX_FILENAME_LENGTH];
   char   filename_output_dir_horizontal_ghosts[MAX_FILENAME_LENGTH];
   char   filename_cat_root[MAX_FILENAME_LENGTH];
   char   filename_groups[MAX_FILENAME_LENGTH];
   char   filename_subgroups[MAX_FILENAME_LENGTH];
   int    i_file =0;
   int    n_files=1;
   sprintf(filename_output_dir_horizontal,       "%s/horizontal",    filename_output_dir);
   sprintf(filename_output_dir_horizontal_ghosts,"%s/ghost_catalogs",filename_output_dir_horizontal);
   mkdir(filename_output_dir,                  02755);
   mkdir(filename_output_dir_horizontal,       02755);
   mkdir(filename_output_dir_horizontal_ghosts,02755);
   strcpy(filename_cat_root,filename_cat);
   strip_path(filename_cat_root);
   sprintf(filename_groups,   "%s/ghosts_%03d.catalog_groups_properties",   filename_output_dir_horizontal_ghosts,j_write);
   sprintf(filename_subgroups,"%s/ghosts_%03d.catalog_subgroups_properties",filename_output_dir_horizontal_ghosts,j_write);
   SID_fopen(filename_groups,"w",&fp_groups);
   SID_fwrite(&i_file,        sizeof(int),1,&fp_groups);
   SID_fwrite(&n_files,       sizeof(int),1,&fp_groups);
   SID_fwrite(&n_group_ghosts,sizeof(int),1,&fp_groups);
   SID_fwrite(&n_group_ghosts,sizeof(int),1,&fp_groups);
   SID_fopen(filename_subgroups,"w",&fp_subgroups);
   SID_fwrite(&i_file,           sizeof(int),1,&fp_subgroups);
   SID_fwrite(&n_files,          sizeof(int),1,&fp_subgroups);
   SID_fwrite(&n_subgroup_ghosts,sizeof(int),1,&fp_subgroups);
   SID_fwrite(&n_subgroup_ghosts,sizeof(int),1,&fp_subgroups);

   // Scan through the trees, writing ghosts in the order that we find them
   int i_group;
   int i_subgroup;
   int flag_null_ghost;
   for(i_group=0,i_subgroup=0;i_group<(n_groups+n_group_ghosts);i_group++){

      // Write ghost groups
      flag_null_ghost=FALSE;
      if(i_group>=n_groups){
         tree_horizontal_ghost_group_info *group_i;
         halo_properties_info *group_initial;
         halo_properties_info *group_final;
         double                time_initial;
         double                time_final;
         double                time_return;
         // Check if this is a null group.  If it is, write it when we write the substructure (below).
         group_i         =&(groups[i_group]);
         flag_null_ghost=check_mode_for_flag(group_i->type,TREE_CASE_GHOST_NULL);
         if(!flag_null_ghost){
            group_initial   =group_properties[(group_i->interp.file_start)%n_wrap][group_i->interp.index_start];
            group_final     =group_properties[(group_i->interp.file_stop)%n_wrap][group_i->interp.index_stop];
            time_initial    =group_i->interp.time_start;
            time_final      =group_i->interp.time_stop;
            time_return     =deltat_a(cosmo,a_list[l_write],0.)/S_PER_YEAR;
            interpolate_halo_local(group_initial,
                                   group_final,
                                   group_return,
                                   time_initial,
                                   time_final,
                                   time_return,
                                   i_group,
                                   group_i->type,
                                   group_i->file_offset,
                                   "");
            SID_fwrite(group_return,sizeof(halo_properties_info),1,&fp_groups);
         }
      }

      // Write the substructure information
      int n_subgroups_i;
      tree_horizontal_ghost_subgroup_info *current;
      current      =groups[i_group].first_substructure;
      n_subgroups_i=0;
      while(current!=NULL){
         if(check_mode_for_flag(current->type,TREE_CASE_GHOST)){
            halo_properties_info *subgroup_initial;
            halo_properties_info *subgroup_final;
            double                time_initial;
            double                time_final;
            double                time_return;
            subgroup_initial=subgroup_properties[(current->interp.file_start)%n_wrap][current->interp.index_start];
            subgroup_final  =subgroup_properties[(current->interp.file_stop)%n_wrap][current->interp.index_stop];
            time_initial    =current->interp.time_start;
            time_final      =current->interp.time_stop;
            time_return     =deltat_a(cosmo,a_list[l_write],0.)/S_PER_YEAR;
            interpolate_halo_local(subgroup_initial,
                                   subgroup_final,
                                   subgroup_return,
                                   time_initial,
                                   time_final,
                                   time_return,
                                   i_group,current->type,
                                   current->file_offset,
                                   "sub");

            // Write the group properties of a null ghost
            n_subgroups_i++;
            if(flag_null_ghost){
               if(n_subgroups_i>1)
                  SID_trap_error("Multiple subgroups have been assigned to a null ghost.",ERROR_LOGIC);
               SID_fwrite(subgroup_return,sizeof(halo_properties_info),1,&fp_groups);
            }

            // Write subgroup ghost properties
            SID_fwrite(subgroup_return,sizeof(halo_properties_info),1,&fp_subgroups);
         }
         current=current->next_substructure;
      }
   }

   // Clean-up
   SID_fclose(&fp_groups);
   SID_fclose(&fp_subgroups);
   SID_free(SID_FARG group_return);
   SID_free(SID_FARG subgroup_return);

   SID_log("Done.",SID_LOG_CLOSE);
}

