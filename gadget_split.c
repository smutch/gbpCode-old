#define  _MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpSPH.h>

int main(int argc, char *argv[]){
  SID_init(&argc,&argv,NULL);

  // Parse command line
  int  snapshot;
  char filename_in_root[MAX_FILENAME_LENGTH];
  char filename_out_root[MAX_FILENAME_LENGTH];
  int  n_split;
  if(argc!=5){
    fprintf(stderr,"\n syntax: %s filename_in snapshot filename_out n_files\n",argv[0]);
    fprintf(stderr," ------\n\n");
    return(ERROR_SYNTAX);
  }
  else{
    strcpy(filename_in_root, argv[1]);
    snapshot=atoi(argv[2]);
    strcpy(filename_out_root,argv[3]);
    n_split=atoi(argv[4]);
  }

  SID_log("Splitting GADGET binary file(s) of {%s} into %d parts...",SID_LOG_OPEN|SID_LOG_TIMER,filename_in_root,n_split);

  // Initialize the plist datastructure and read the header
  plist_info plist;
  init_plist(&plist,NULL,GADGET_LENGTH,GADGET_MASS,GADGET_VELOCITY);

  // Read the header and determine the input file-format
  gadget_header_info header;
  int                i_file;
  int                flag_filefound;
  int                flag_multifile;
  int                flag_file_type;
  flag_filefound=init_gadget_read(filename_in_root,snapshot,&flag_multifile,&flag_file_type,&header);
  if(!flag_filefound)
     SID_trap_error("File not found.",ERROR_LOGIC);

  // Loop over the number of files in the snapshot
  int n_files;
  if(flag_filefound)
     n_files=header.n_files;
  else
     n_files=0;
  if(n_files<SID.n_proc)
     SID_log("You are processing %d files with %d cores.  %d cores will go unused.",SID_LOG_COMMENT,n_files,SID.n_proc,SID.n_proc-n_files);
  char filename_in[MAX_FILENAME_LENGTH];
  char filename_out[MAX_FILENAME_LENGTH];
  for(i_file=0;i_file<n_files;i_file++){
     if(n_files>1)
        SID_log("Processing file(s) %d->%d...",SID_LOG_OPEN,i_file,MIN(i_file+SID.n_proc-1,n_files-1));

     // Create the directory needed for the snapshot
     set_gadget_filename(filename_out_root,snapshot,0,TRUE,0,filename_out);
     if(i_file==0){
        FILE *fp_test;
        char  filename_out_directory[MAX_FILENAME_LENGTH];
        strcpy(filename_out_directory,filename_out);
        strip_file_root(filename_out_directory);
        if((fp_test=fopen(filename_out_directory,"r"))==NULL)
           mkdir(filename_out_directory,02755);
        else
           fclose(fp_test);
     }

     // Open input file and read header
     int   n_in[N_GADGET_TYPE];
     int   n_in_total=0;
     int   n_written[N_GADGET_TYPE];
     FILE *fp_in;
     FILE *fp_out;
     int   record_length_in;
     int   record_length_out;
     int   record_length_header;
     int   i_type;
     set_gadget_filename(filename_in_root,snapshot,i_file,flag_multifile,flag_file_type,filename_in);
     if((fp_in=fopen(filename_in,"r"))==NULL)
        SID_trap_error("Could not open {%s}.",ERROR_IO_OPEN,filename_in);
     fread(&record_length_in,sizeof(int),1,fp_in);
     fread(&header,sizeof(gadget_header_info),1,fp_in);
     fread(&record_length_out,sizeof(int),1,fp_in);
     for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
       n_in[i_type]      =header.n_file[i_type];
       n_in_total       +=header.n_file[i_type];
       n_written[i_type] =0;
       header.n_file[i_type]=0;
     }
     header.n_files*=n_split;
     record_length_header=record_length_out;

     int i_split;
     int j_split;
     for(i_split=i_file*n_split,j_split=0;j_split<n_split;i_split++,j_split++){
        // Rewind to the beginning of the input file and then skip the header
        rewind(fp_in);
        fseeko(fp_in,4+sizeof(gadget_header_info)+4,SEEK_SET);

        // Open new file for writing
        char filename_in[MAX_FILENAME_LENGTH];
        set_gadget_filename(filename_out_root,snapshot,i_split,TRUE,0,filename_out);
        fp_out=fopen(filename_out,"w");

        // Set the particle counts for this file
        int n_particles_i=0;
        int n_buffer_max =0;
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
           if(j_split==(n_split-1))
              header.n_file[i_type]=n_in[i_type]-n_written[i_type];
           else
              header.n_file[i_type]=n_in[i_type]/n_split;
           n_particles_i+=header.n_file[i_type];
           n_buffer_max  =MAX(n_buffer_max,header.n_file[i_type]);
        }

        // Create read buffer (this is long enough for all blocks, even if IDs are long)
        char *buffer;
        buffer=(char *)SID_malloc(3*sizeof(GBPREAL)*n_buffer_max);

        // Write the header
        fwrite(&record_length_header,sizeof(int),1,fp_out);
        fwrite(&header,sizeof(gadget_header_info),1,fp_out);
        fwrite(&record_length_header,sizeof(int),1,fp_out);

        // Read/Write Positions
        int record_length_write;
        record_length_write=3*sizeof(GBPREAL)*n_particles_i;
        fwrite(&record_length_write,sizeof(int),1,fp_out);
        fread(&record_length_in,    sizeof(int),1,fp_in);
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
          fseeko(fp_in,3*sizeof(GBPREAL)*n_written[i_type],SEEK_CUR);
          fread(buffer,3*header.n_file[i_type],sizeof(GBPREAL),fp_in);
          fseeko(fp_in,3*sizeof(GBPREAL)*(n_in[i_type]-n_written[i_type]-header.n_file[i_type]),SEEK_CUR);
          fwrite(buffer,3*sizeof(GBPREAL)*header.n_file[i_type],1,fp_out);
        }
        fread(&record_length_out,sizeof(int),1,fp_in);
        fwrite(&record_length_write,sizeof(int),1,fp_out);

        // Read/Write Velocities
        fwrite(&record_length_write,sizeof(int),1,fp_out);
        fread(&record_length_in,    sizeof(int),1,fp_in);
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
          fseeko(fp_in,3*sizeof(GBPREAL)*n_written[i_type],SEEK_CUR);
          fread(buffer,3*header.n_file[i_type],sizeof(GBPREAL),fp_in);
          fseeko(fp_in,3*sizeof(GBPREAL)*(n_in[i_type]-n_written[i_type]-header.n_file[i_type]),SEEK_CUR);
          fwrite(buffer,3*sizeof(GBPREAL)*header.n_file[i_type],1,fp_out);
        }
        fread(&record_length_out,sizeof(int),1,fp_in);
        fwrite(&record_length_write,sizeof(int),1,fp_out);

        // Read the IDs header length and use it's value to determine if IDs are long or short
        int id_byte_size;
        fread(&record_length_in,sizeof(int),1,fp_in);
        if(record_length_in/sizeof(int)==n_in_total)
           id_byte_size=sizeof(int);
        else if(record_length_in/sizeof(long long)==n_in_total)
           id_byte_size=sizeof(long long);
        else
           SID_trap_error("Invalid input IDs block size (%d).",ERROR_LOGIC,record_length_in);

        // Read/Write IDs
        record_length_write=n_particles_i*id_byte_size;
        fwrite(&record_length_write,sizeof(int),1,fp_out);
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++){
           fseeko(fp_in,id_byte_size*n_written[i_type],SEEK_CUR);
           fread(buffer,header.n_file[i_type],id_byte_size,fp_in);
           fseeko(fp_in,id_byte_size*(n_in[i_type]-n_written[i_type]-header.n_file[i_type]),SEEK_CUR);
           fwrite(buffer,id_byte_size*header.n_file[i_type],1,fp_out);
        }
        fread(&record_length_out,   sizeof(int),1,fp_in);
        fwrite(&record_length_write,sizeof(int),1,fp_out);

        // Update the particle counters
        for(i_type=0;i_type<N_GADGET_TYPE;i_type++)
           n_written[i_type]+=header.n_file[i_type];

        // Clean-up
        SID_free(SID_FARG buffer);
     }
     fclose(fp_in);
     fclose(fp_out);
     if(n_files>1)
        SID_log("Done.",SID_LOG_CLOSE);
  }
  
  // Clean-up 
  free_plist(&plist);
  SID_log("Done.",SID_LOG_CLOSE);

  SID_exit(ERROR_NONE);
}

