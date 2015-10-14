#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpSPH.h>

void read_smooth(plist_info *plist,
                 char       *filename_root_in,
                 int         snapshot_number,
                 int         mode){
  char    filename[256];
  size_t  n_particles_all_mem;
  size_t  n_particles_local;
  size_t  n_particles_total;
  int     i_quantity;
  int     n_quantities=3;
  int    *used;
  char   *species_name;
  char    var_name[256];
  char    unit_name[256];
  double  unit_factor;
  int     i_file;
  size_t *ids;
  size_t *ids_index;
  int     n_particles_file;
  int     offset;
  int     n_files;
  void      *id_buf;
  size_t    *id_buf_index=NULL;
  int       *id_buf_i;
  long long *id_buf_L;
  long long *mark;
  size_t  i_particle;
  size_t  j_particle;
  void   *buffer;
  float  *local_array;
  int     n_mark;
  float  *r_smooth_array;  
  float  *rho_array;  
  float  *sigma_v_array;  
  char   *read_array;  
  double  expansion_factor;
  double  h_Hubble;
  int     flag_filefound=FALSE;
  int     flag_multifile=FALSE;
  int     flag_file_type;
  int     flag_LONGIDs;
  int     read_rank=MASTER_RANK;
  int     flag_log_sigma=FALSE;
  int     flag_log_rho=FALSE;
  FILE   *fp;

  SID_log("Reading smooth file {%s}...",SID_LOG_OPEN|SID_LOG_TIMER,filename_root_in);

  // Read header info
  SID_log("Reading header information...",SID_LOG_OPEN);
  smooth_header_info header;
  flag_filefound   =init_smooth_read(filename_root_in,snapshot_number,&flag_multifile,&flag_file_type,&header);
  SID_log("n_files    =%d",  SID_LOG_COMMENT,header.n_files);
  SID_log("n_particles=%lld",SID_LOG_COMMENT,header.n_particles_total);
  SID_log("Done.",SID_LOG_CLOSE);

  // Interpret the mode passed to this function
  int flag_logs_used   =FALSE;
  int flag_log_quantity=FALSE;
  if(check_mode_for_flag(mode,READ_SMOOTH_LOG_SIGMA)){
     flag_log_sigma=TRUE;
     flag_logs_used=TRUE;
  }
  else
     flag_log_sigma=FALSE;
  if(check_mode_for_flag(mode,READ_SMOOTH_LOG_RHO)){
     flag_log_rho  =TRUE;
     flag_logs_used=TRUE;
  }
  else
     flag_log_rho=FALSE;

  // A file was found ... 
  if(flag_filefound){
     // Communicate the header
     n_particles_file =header.n_particles_file;
     offset           =header.offset;
     n_particles_total=header.n_particles_total;
     n_files          =MAX(1,header.n_files);
     SID_Bcast(&n_particles_file, (int)sizeof(int),      read_rank,SID.COMM_WORLD);
     SID_Bcast(&offset,           (int)sizeof(int),      read_rank,SID.COMM_WORLD);
     SID_Bcast(&n_particles_total,(int)sizeof(long long),read_rank,SID.COMM_WORLD);
     SID_Bcast(&n_files,          (int)sizeof(int),      read_rank,SID.COMM_WORLD);

     // Fetch the number of particles and their ids
     species_name     =plist->species[GADGET_TYPE_DARK];
     n_particles_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",    species_name))[0]; 
     n_particles_all_mem  =((size_t *)ADaPS_fetch(plist->data,"n_all_%s",species_name))[0]; 
     ids              = (size_t *)ADaPS_fetch(plist->data,"id_%s",   species_name); 
     expansion_factor =((double *)ADaPS_fetch(plist->data,"expansion_factor"))[0]; 
     h_Hubble         =((double *)ADaPS_fetch(plist->data,"h_Hubble"))[0]; 

     // Check if we are dealing with LONG IDs or not
     if(ADaPS_exist(plist->data,"flag_LONGIDs"))
       flag_LONGIDs=TRUE;
     else
       flag_LONGIDs=FALSE;

     // Sort particle IDs
     SID_log("Sorting particle IDs...",SID_LOG_OPEN|SID_LOG_TIMER);
     merge_sort(ids,(size_t)n_particles_local,&ids_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
     SID_log("Done.",SID_LOG_CLOSE);

     // Allocate arrays
     SID_log("Allocating arrays...",SID_LOG_OPEN);
     read_array=(char *)SID_calloc(sizeof(char)*n_particles_local);
     for(i_quantity=0;i_quantity<n_quantities;i_quantity++){
       switch(i_quantity){
       case 0:
         r_smooth_array=(float *)SID_malloc(sizeof(float)*n_particles_local);
         ADaPS_store(&(plist->data),r_smooth_array,"r_smooth_%s",ADaPS_DEFAULT,species_name);
         break;
       case 1:
         rho_array=(float *)SID_malloc(sizeof(float)*n_particles_local);
         ADaPS_store(&(plist->data),rho_array,"rho_%s",ADaPS_DEFAULT,species_name);
         break;
       case 2:
         sigma_v_array=(float *)SID_malloc(sizeof(float)*n_particles_local);
         ADaPS_store(&(plist->data),sigma_v_array,"sigma_v_%s",ADaPS_DEFAULT,species_name);
         break;
       }
     }
     SID_log("Done.",SID_LOG_CLOSE);

     // Read each file in turn
     pcounter_info pcounter;
     size_t n_particles_read=0;
     SID_log("Performing read...",SID_LOG_OPEN|SID_LOG_TIMER);
     SID_init_pcounter(&pcounter,n_particles_total,10);
     for(i_file=0;i_file<n_files;i_file++){
       set_smooth_filename(filename_root_in,snapshot_number,i_file,flag_multifile,flag_file_type,filename);
       if((fp=fopen(filename,"r"))!=NULL){
          fread_verify(&(header.n_particles_file), sizeof(int),      1,fp);
          fread_verify(&(header.offset),           sizeof(int),      1,fp);
          fread_verify(&(header.n_particles_total),sizeof(long long),1,fp);
          fread_verify(&(header.n_files),          sizeof(int),      1,fp);
       }
       else
          SID_trap_error("Could not open smooth file {%s}",ERROR_IO_OPEN,filename);
       n_particles_file =header.n_particles_file;
       offset           =header.offset;
       n_particles_total=header.n_particles_total;
       n_files          =MAX(1,header.n_files);
       SID_Bcast(&n_particles_file, (int)sizeof(int),      read_rank,SID.COMM_WORLD);
       SID_Bcast(&offset,           (int)sizeof(int),      read_rank,SID.COMM_WORLD);
       SID_Bcast(&n_particles_total,(int)sizeof(long long),read_rank,SID.COMM_WORLD);
       SID_Bcast(&n_files,          (int)sizeof(int),      read_rank,SID.COMM_WORLD);

       // Read IDs
       if(flag_LONGIDs){
         if(i_file==0) SID_log("(long long) IDs...",SID_LOG_CONTINUE);
         id_buf  =SID_malloc(sizeof(long long)*n_particles_file);
         id_buf_L=(long long *)id_buf;
         if(SID.My_rank==read_rank){
           fseeko(fp,(size_t)(3*n_particles_file*sizeof(float)),SEEK_CUR);
           fread_verify(id_buf,sizeof(long long),n_particles_file,fp);
         }
         SID_Barrier(SID.COMM_WORLD);
         SID_Bcast(id_buf_L,(int)(n_particles_file*sizeof(long long)),read_rank,SID.COMM_WORLD);
         merge_sort(id_buf_L,(size_t)n_particles_file,&id_buf_index,SID_SIZE_T,SORT_COMPUTE_INDEX,FALSE);
       }
       else{
         if(i_file==0) SID_log("(int) IDs...",SID_LOG_CONTINUE);
         id_buf  =SID_malloc(sizeof(int)*n_particles_file);
         id_buf_i=(int *)id_buf;
         if(SID.My_rank==read_rank){
           fseeko(fp,(size_t)(3*n_particles_file*sizeof(float)),SEEK_CUR);
           fread_verify(id_buf,sizeof(int),n_particles_file,fp);
         }
         SID_Barrier(SID.COMM_WORLD);
         SID_Bcast(id_buf_i,(int)(n_particles_file*sizeof(int)),read_rank,SID.COMM_WORLD);
         merge_sort(id_buf_i,(size_t)n_particles_file,&id_buf_index,SID_INT,SORT_COMPUTE_INDEX,FALSE);
       }

       // Create local particle mapping
       int n_mark_i=0;
       mark=(long long *)SID_malloc(sizeof(long long)*n_particles_file);
       for(i_particle=0;i_particle<n_particles_file;i_particle++) 
         mark[i_particle]=-1;
       if(flag_LONGIDs){
         for(i_particle=0,j_particle=0,n_mark=0;i_particle<n_particles_file && j_particle<n_particles_local;i_particle++){
           while(j_particle<n_particles_local-1 && ids[ids_index[j_particle]]<id_buf_L[id_buf_index[i_particle]]) j_particle++;
           if(ids[ids_index[j_particle]]==id_buf_L[id_buf_index[i_particle]]){
             mark[id_buf_index[i_particle]]=(long long)ids_index[j_particle];
             n_mark++;
             n_mark_i++;
           }
         }
       }
       else{
         for(i_particle=0,j_particle=0,n_mark=0;i_particle<n_particles_file && j_particle<n_particles_local;i_particle++){
           while(j_particle<n_particles_local-1 && ids[ids_index[j_particle]]<id_buf_i[id_buf_index[i_particle]]) j_particle++;
           if(ids[ids_index[j_particle]]==id_buf_i[id_buf_index[i_particle]]){
             mark[id_buf_index[i_particle]]=(long long)ids_index[j_particle];
             n_mark++;
             n_mark_i++;
           }
         }
       }
       SID_free(SID_FARG id_buf);
       SID_free(SID_FARG id_buf_index);

       // Move to the start of the particle quantities
       if(SID.My_rank==read_rank){
         rewind(fp);
         fread_verify(&n_particles_file, sizeof(int),      1,fp);
         fread_verify(&offset,           sizeof(int),      1,fp);
         fread_verify(&n_particles_total,sizeof(long long),1,fp);
         fread_verify(&n_files,          sizeof(int),      1,fp);
         n_files=MAX(1,n_files);
       }
       SID_Bcast(&n_particles_file, (int)sizeof(int),      read_rank,SID.COMM_WORLD);
       SID_Bcast(&offset,           (int)sizeof(int),      read_rank,SID.COMM_WORLD);
       SID_Bcast(&n_particles_total,(int)sizeof(long long),read_rank,SID.COMM_WORLD);
       SID_Bcast(&n_files,          (int)sizeof(int),      read_rank,SID.COMM_WORLD);
       buffer=SID_malloc(sizeof(float)*n_particles_file);
       for(i_quantity=0;i_quantity<n_quantities;i_quantity++){
         int flag_log_quantity;
         switch(i_quantity){
         case 0:
           sprintf(var_name,"r_smooth_%s",species_name);
           sprintf(unit_name,"Mpc");
           unit_factor=plist->length_unit/h_Hubble;
           local_array=r_smooth_array;
           break;
         case 1:
           sprintf(var_name,"rho_%s",species_name);
           sprintf(unit_name,"Msol/Mpc^3");
           unit_factor=h_Hubble*h_Hubble*plist->mass_unit/pow(plist->length_unit,3.);
           local_array=rho_array;
           break;
         case 2:
           sprintf(var_name,"sigma_v_%s",species_name);
           sprintf(unit_name,"km/s");
           unit_factor=sqrt(expansion_factor)*plist->velocity_unit;
           local_array=sigma_v_array;
           break;
         }

         // Read next quantity
         if(SID.My_rank==read_rank)
           fread_verify(buffer,sizeof(float),n_particles_file,fp);
         SID_Barrier(SID.COMM_WORLD);
         SID_Bcast(buffer,(int)(n_particles_file*sizeof(float)),read_rank,SID.COMM_WORLD);

         // Place in final array
         if(n_mark_i>0){
           if(i_quantity==0){
             for(i_particle=0;i_particle<n_particles_file;i_particle++)
               if(mark[i_particle]>=0) read_array[mark[i_particle]]=TRUE;
           }
           for(i_particle=0;i_particle<n_particles_file;i_particle++)
             if(mark[i_particle]>=0) local_array[mark[i_particle]]=((float *)buffer)[i_particle]*unit_factor;
         }
       }
       SID_free(SID_FARG mark);
       SID_free(SID_FARG buffer);
       if(SID.My_rank==read_rank)
         fclose(fp);
       n_particles_read+=n_particles_file;
       if(n_files>1)
          SID_check_pcounter(&pcounter,n_particles_read);
     } // i_file
     SID_free(SID_FARG ids_index);
     SID_Barrier(SID.COMM_WORLD);
     SID_log("Done.",SID_LOG_CLOSE);

     // Check that all particles have been treated
     size_t sum_check=0;
     for(i_particle=0;i_particle<n_particles_local;i_particle++)
        sum_check+=(size_t)read_array[i_particle];
     SID_Allreduce(SID_IN_PLACE,&sum_check,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
     if(sum_check!=n_particles_all_mem)
        SID_trap_error("Only %lld of %lld particles were set.",ERROR_LOGIC,sum_check,n_particles_all_mem);
     SID_free(SID_FARG read_array);

     SID_log("Summary...",SID_LOG_OPEN);
     for(i_quantity=0;i_quantity<n_quantities;i_quantity++){
       char var_name[32];
       switch(i_quantity){
       case 0:
         sprintf(var_name,"Lengths:   ");
         sprintf(unit_name,"Mpc");
         unit_factor=1./M_PER_MPC;
         local_array=r_smooth_array;
         break;
       case 1:
         sprintf(var_name,"Densities: ");
         sprintf(unit_name,"Msol/Mpc^3");
         unit_factor=M_PER_MPC*M_PER_MPC*M_PER_MPC/M_SOL;
         local_array=rho_array;
         break;
       case 2:
         sprintf(var_name,"sigmas_v's:");
         sprintf(unit_name,"km/s");
         unit_factor=1e-3;
         local_array=sigma_v_array;
         break;
       }

       // Report some stats
       float var_min,var_max,var_mean;
       calc_min_global(local_array,
                       &var_min,
                       n_particles_local,
                       SID_FLOAT,
                       CALC_MODE_DEFAULT,
                       SID.COMM_WORLD);
       calc_max_global(local_array,
                       &var_max,
                       n_particles_local,
                       SID_FLOAT,
                       CALC_MODE_DEFAULT,
                       SID.COMM_WORLD);
       calc_mean_global(local_array,
                        &var_mean,
                        n_particles_local,
                        SID_FLOAT,
                        CALC_MODE_DEFAULT,
                        SID.COMM_WORLD);

       // Remove NaN's from sigma_v's if needed
       size_t n_NaN=0;
       for(i_particle=0;i_particle<n_particles_local;i_particle++){
          if(isnan(local_array[i_particle])){
             local_array[i_particle]=1000.*0.5*(var_min+var_max);
             n_NaN++;
          }
       }

       if(n_NaN>0)
          SID_log("%s min=%e max=%e mean=%e [%s] (%lld are NaN)",
                  SID_LOG_COMMENT,
                  var_name,
                  var_min*unit_factor,
                  var_max*unit_factor,
                  var_mean*unit_factor,
                  unit_name,
                  n_NaN);
       else
          SID_log("%s min=%e max=%e mean=%e [%s]",
                  SID_LOG_COMMENT,
                  var_name,
                  var_min*unit_factor,
                  var_max*unit_factor,
                  var_mean*unit_factor,
                  unit_name);
     }
     SID_log("",SID_LOG_CLOSE|SID_LOG_NOPRINT);

     if(flag_logs_used){
        SID_log("Taking needed logs of quantities...",SID_LOG_OPEN|SID_LOG_TIMER);   
        for(i_quantity=0;i_quantity<n_quantities;i_quantity++){
          switch(i_quantity){
          case 0:
            unit_factor=1./M_PER_MPC;
            local_array=r_smooth_array;
            flag_log_quantity=FALSE;
            break;
          case 1:
            unit_factor=M_PER_MPC*M_PER_MPC*M_PER_MPC/M_SOL;
            local_array=rho_array;
            flag_log_quantity=flag_log_rho;
            break;
          case 2:
            unit_factor=1e-3;
            local_array=sigma_v_array;
            flag_log_quantity=flag_log_sigma;
            break;
          }
          if(flag_log_quantity){
             for(i_particle=0;i_particle<n_particles_local;i_particle++)
                local_array[i_particle]=take_log10(local_array[i_particle]);
          }
        }
        SID_log("Done.",SID_LOG_CLOSE);
     }
     SID_Barrier(SID.COMM_WORLD);
     SID_log("Done.",SID_LOG_CLOSE);
  }
  else
    SID_trap_error("Could not find file with root {%s}",ERROR_IO_OPEN,filename_root_in);
}

