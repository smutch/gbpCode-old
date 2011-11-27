#define _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  plist_info plist;
  char       filename_in_root[256];
  char       filename_out_root[256];
  char       filename_in[256];
  char       filename_out[256];
  int        subsample_factor;
  int        i_snap,i_file,j_file;
  int        n_files;
  unsigned int n_local[N_GADGET_TYPE];
  int          i_type;
  size_t     i_particle,j_particle;
  FILE      *fp;
  FILE      *fp_out;
  int        record_length_in;
  int        record_length_out;
  GBPREAL       position[3];
  GBPREAL       velocity[3];
  int        ID_int;
  RNG_info   RNG;
  int        seed_init=10463743;
  int        seed;
  GBPREAL       RNG_max;
  gadget_header_info header;
  long long n_all_keep[N_GADGET_TYPE];
  long long n_all_init[N_GADGET_TYPE];
  unsigned int n_all_high_word[N_GADGET_TYPE];
  unsigned int n_all[N_GADGET_TYPE];
  double    mass_array[N_GADGET_TYPE];
  int       n_keep_total;
  int       n_init_total;
  int       record_length;
  int       flag_long_ids;
  int       id_in_int;
  long long id_in_long;
  int       n_groups;
  int       flag_filefound;
  int       flag_multifile;
  int       flag_file_type;

  SID_init(&argc,&argv,NULL);

  // Parse command line
  if(argc!=4){
    fprintf(stderr,"\n syntax: %s filename_in_root snapshot_number subsample_factor\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in_root, argv[1]);
    i_snap          =atoi(argv[2]);
    subsample_factor=atoi(argv[3]);
  }
  sprintf(filename_out_root,"%s_subsample",filename_in_root);

  RNG_max=1./(GBPREAL)(subsample_factor);

  SID_log("Subsampling GADGET file by a factor of %d...",SID_LOG_OPEN,subsample_factor);

  // Grab info from first header
  flag_filefound=init_gadget_read(filename_in_root,i_snap,&flag_multifile,&flag_file_type,&header);
  if(flag_filefound)
    n_files =header.n_files;
  else
    n_files=0;
  if(n_files>0){
  init_seed_from_clock(&seed_init);
  SID_log("Seed=%d",SID_LOG_COMMENT,seed_init);
  
  SID_log("Counting particles in %d files...",SID_LOG_OPEN,n_files);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
    n_all_keep[i_type]=0;
  for(i_file=SID.My_rank,j_file=0;j_file<n_files;i_file+=SID.n_proc,j_file+=SID.n_proc){
    SID_log("Processing file(s) %d->%d...",SID_LOG_OPEN,j_file,MIN(j_file+SID.n_proc-1,n_files-1));
    if(i_file<n_files){
      set_gadget_filename(filename_in_root,i_snap,i_file,flag_multifile,flag_file_type,filename_in);
      fp=fopen(filename_in,"r");
      seed=seed_init+1387*i_file;
      init_RNG(&seed,&RNG,RNG_DEFAULT);
      for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
        for(j_particle=0;j_particle<header.n_file[i_type];j_particle++){
          if(random_number(&RNG)<=RNG_max)
            n_all_keep[i_type]++;
        }
      }
      free_RNG(&RNG);
      fclose(fp);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  SID_Allreduce(SID_IN_PLACE,n_all_keep,N_GADGET_TYPE,SID_LONG_LONG,SID_SUM,SID.COMM_WORLD);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_log("Header properties will be...",SID_LOG_OPEN);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
    n_all_high_word[i_type]=(unsigned int)(((long long)n_all_keep[i_type])>>32);
    n_all[i_type]          =(unsigned int)(((long long)n_all_keep[i_type])-(((long long)(n_all_high_word[i_type]))<<32));
    n_all_init[i_type]     =((long long)(header.n_all[i_type]))+(((long long)(header.n_all_high_word[i_type]))<<32);
    if(n_all_keep[i_type]>0)
      mass_array[i_type]=header.mass_array[i_type]*((double)n_all_init[i_type]/(double)n_all_keep[i_type]);
    else
      mass_array[i_type]=0.;
  }
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
    SID_log("mass_array[%d]=%le (was %le)", SID_LOG_COMMENT,i_type,mass_array[i_type],header.mass_array[i_type]);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
    SID_log("n_all[%d]     =%d (was %d)",SID_LOG_COMMENT,i_type,n_all[i_type],header.n_all[i_type]);
  for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
    SID_log("high_word[%d] =%d (was %d)",SID_LOG_COMMENT,i_type,n_all_high_word[i_type],header.n_all_high_word[i_type]);
  SID_log("Done.",SID_LOG_CLOSE);
  for(i_file=SID.My_rank,j_file=0;j_file<n_files;i_file+=SID.n_proc,j_file+=SID.n_proc){
    SID_log("Processing file(s) %d->%d...",SID_LOG_OPEN,j_file,MIN(j_file+SID.n_proc-1,n_files-1));
    if(i_file<n_files){
      set_gadget_filename(filename_in_root, i_snap,i_file,flag_multifile,flag_file_type,filename_in);
      set_gadget_filename(filename_out_root,i_snap,i_file,flag_multifile,flag_file_type,filename_out);
      fp    =fopen(filename_in, "r");
      fp_out=fopen(filename_out,"w");

      // Header
      fread(&record_length_in,sizeof(int),1,fp);
      fread(&header,sizeof(gadget_header_info),1,fp);
      fread(&record_length_out,sizeof(int),1,fp);
      for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
        n_local[i_type]               =header.n_file[i_type];
        header.n_file[i_type]         =0;
        header.mass_array[i_type]     =mass_array[i_type];
        header.n_all[i_type]          =n_all[i_type];
        header.n_all_high_word[i_type]=n_all_high_word[i_type];
      }
      seed=seed_init+1387*i_file;
      init_RNG(&seed,&RNG,RNG_DEFAULT);
      for(i_type=0,n_keep_total=0,n_init_total=0;i_type<N_GADGET_TYPE;i_type++){
        for(j_particle=0;j_particle<n_local[i_type];j_particle++,i_particle++){
          if(random_number(&RNG)<=RNG_max){
            header.n_file[i_type]++;
            n_keep_total++;
          }
          n_init_total++;
        }
      }
      free_RNG(&RNG);
      record_length=record_length_in;
      fwrite(&record_length,sizeof(int),1,fp_out);
      fwrite(&header,sizeof(gadget_header_info),1,fp_out);
      fwrite(&record_length,sizeof(int),1,fp_out);

      // Positions
      record_length=3*sizeof(GBPREAL)*n_keep_total;
      seed=seed_init+1387*i_file;
      init_RNG(&seed,&RNG,RNG_DEFAULT);
      fread(&record_length_in,sizeof(int),1,fp);
      fwrite(&record_length,sizeof(int),1,fp_out);
      for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
        for(j_particle=0;j_particle<n_local[i_type];j_particle++){
          fread(position,sizeof(GBPREAL),3,fp);
          if(random_number(&RNG)<=RNG_max)
            fwrite(position,sizeof(GBPREAL),3,fp_out);
        }
      }
      fread(&record_length_out,sizeof(int),1,fp);
      fwrite(&record_length,sizeof(int),1,fp_out);
      free_RNG(&RNG);

      // Velocities
      record_length=3*sizeof(GBPREAL)*n_keep_total;
      seed=seed_init+1387*i_file;
      init_RNG(&seed,&RNG,RNG_DEFAULT);
      fread(&record_length_in,sizeof(int),1,fp);
      fwrite(&record_length,sizeof(int),1,fp_out);
      for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
        for(j_particle=0;j_particle<n_local[i_type];j_particle++){
          fread(velocity,sizeof(GBPREAL),3,fp);
          if(random_number(&RNG)<=RNG_max)
            fwrite(velocity,sizeof(GBPREAL),3,fp_out);
        }
      }
      fread(&record_length_out,sizeof(int),1,fp);
      fwrite(&record_length,sizeof(int),1,fp_out);
      free_RNG(&RNG);

      // IDs
      seed=seed_init+1387*i_file;
      init_RNG(&seed,&RNG,RNG_DEFAULT);
      fread(&record_length_in,sizeof(int),1,fp);
      if(record_length_in/sizeof(int)==n_init_total){
        flag_long_ids=FALSE;
        record_length=n_keep_total*sizeof(int);
      }
      else{
        flag_long_ids=TRUE;
        record_length=n_keep_total*sizeof(long long);
      }
      fwrite(&record_length,sizeof(int),1,fp_out);
      if(flag_long_ids){
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
          for(j_particle=0;j_particle<n_local[i_type];j_particle++){
            fread(&id_in_long,sizeof(long long),1,fp);
            if(random_number(&RNG)<=RNG_max)
              fwrite(&id_in_long,sizeof(long long),1,fp_out);
          }
        }
      }
      else{
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
          for(j_particle=0;j_particle<n_local[i_type];j_particle++){
            fread(&id_in_int,sizeof(int),1,fp);
            if(random_number(&RNG)<=RNG_max)
              fwrite(&id_in_int,sizeof(int),1,fp_out);
          }
        }
      }
      fread(&record_length_out,sizeof(int),1,fp);
      fwrite(&record_length,sizeof(int),1,fp_out);
      free_RNG(&RNG);

      fclose(fp);
      fclose(fp_out);
    }
    SID_log("Done.",SID_LOG_CLOSE);
  }
  }
  else
    SID_log_error("File not found.",ERROR_IO_READ);
  SID_log("Done.",SID_LOG_CLOSE);
  SID_exit(ERROR_NONE);
}
