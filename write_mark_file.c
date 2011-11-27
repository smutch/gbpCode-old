#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <gbpLib.h>
#include <gbpSPH.h>

/* In the future ... receive a list of variables to make this adaptive */

void write_mark_file(plist_info *plist,
                     char       *mark_name,
                     char       *filename_out){
  int     i_species;
  size_t  i_particle;
  int     i_rank;
  size_t  i_mark;
  int    *mark_local;
  size_t *ids_local;
  size_t *mark_list_local;
  size_t  n_mark_total;
  size_t  n_mark_type_local[N_GADGET_TYPE];
  size_t  n_mark_local;
  size_t  n_particle_local;
  SID_fp  fp_mark_file;
  size_t  i_start_local;
  size_t  n_mark_bcast;
  int     n_chunks=10;
  markfile_header_info header={N_GADGET_TYPE};

  // Count the number of particles that we are going to write  
  for(i_species=0,n_mark_local=0;i_species<N_GADGET_TYPE;i_species++){
    n_mark_type_local[i_species]=0;
    if(ADaPS_exist(plist->data,"n_%s",plist->species[i_species])){
      n_particle_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",plist->species[i_species]))[0];
      mark_local      =(int *)ADaPS_fetch(plist->data,"%s_%s",mark_name,plist->species[i_species]);
      for(i_particle=0;i_particle<n_particle_local;i_particle++){
        if(mark_local[i_particle])
          n_mark_type_local[i_species]++;
      }
      n_mark_local+=n_mark_type_local[i_species];
      SID_Allreduce(&(n_mark_type_local[i_species]),&(header.n_mark_species[i_species]),1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);
    }
  }
  SID_Allreduce(&n_mark_local,&n_mark_total,1,SID_SIZE_T,SID_SUM,SID.COMM_WORLD);

  // Open the file and write the header
  SID_fopen_chunked(filename_out,
                    "w",
                    &fp_mark_file,
                    &header,
                    sizeof(markfile_header_info),
                    n_mark_total,
                    sizeof(size_t),
                    n_chunks);

  // Write the particles for each species in turn
  for(i_species=0,n_mark_local=0;i_species<N_GADGET_TYPE;i_species++){
    if(header.n_mark_species[i_species]>0){
       n_particle_local=((size_t *)ADaPS_fetch(plist->data,"n_%s",           plist->species[i_species]))[0];
       mark_local      = (int    *)ADaPS_fetch(plist->data,"%s_%s",mark_name,plist->species[i_species]);
       ids_local       = (size_t *)ADaPS_fetch(plist->data,"id_%s",          plist->species[i_species]);
       if(n_mark_type_local[i_species]>0){
         mark_list_local =(size_t *)SID_calloc(sizeof(size_t)*n_mark_type_local[i_species]);
         for(i_particle=0,i_mark=0;i_particle<n_particle_local;i_particle++){
           if(mark_local[i_particle])
             mark_list_local[i_mark++]=ids_local[i_particle];
         }
       }
       i_start_local=0;
#if USE_MPI
       for(i_rank=0;i_rank<SID.n_proc;i_rank++){
         if(SID.My_rank==i_rank)
           n_mark_bcast=n_mark_type_local[i_species];
         SID_Bcast(&n_mark_bcast,sizeof(size_t),i_rank,SID.COMM_WORLD);
         if(i_rank<SID.My_rank)
           i_start_local+=n_mark_bcast;
       }
#endif
       SID_log("Writing %lld %s particles...",SID_LOG_OPEN,header.n_mark_species[i_species],plist->species[i_species]);
       SID_fwrite_chunked(mark_list_local,
                          n_mark_type_local[i_species],
                          i_start_local,
                          &fp_mark_file);
       SID_log("Done.",SID_LOG_CLOSE);
       SID_free(SID_FARG mark_list_local);
    }
  }
  SID_fclose_chunked(&fp_mark_file);

}
