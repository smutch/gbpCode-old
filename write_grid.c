#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpCosmo.h>
#include <gbpSPH.h>
#include <gbpHalos.h>
#include <gbpClustering.h>

void write_grid(field_info *field,char *filename_out_root,int i_grid,int n_grids,int mass_assignment_scheme,char *grid_identifier,double box_size){
   // Now that all 4 runs are done, let's write the results
   SID_log("Writing {%s} grid...",SID_LOG_OPEN,grid_identifier);

   // Set output filename
   char filename_out[MAX_FILENAME_LENGTH];
   sprintf(filename_out,"%s_grid.dat",filename_out_root);

   // Write header if this is the first grid
   FILE *fp_out=NULL;
   if(SID.I_am_Master){
      if(i_grid==0){
         fp_out=fopen(filename_out,"w");
         fwrite(&(field->n[0]),sizeof(int),   1,fp_out);
         fwrite(&(field->n[1]),sizeof(int),   1,fp_out);
         fwrite(&(field->n[2]),sizeof(int),   1,fp_out);
         fwrite(&(box_size),   sizeof(double),1,fp_out);
         fwrite(&(box_size),   sizeof(double),1,fp_out);
         fwrite(&(box_size),   sizeof(double),1,fp_out);
         fwrite(&(n_grids),    sizeof(int),   1,fp_out);
         switch(mass_assignment_scheme){
            case MAP2GRID_DIST_DWT20:
            case MAP2GRID_DIST_DWT12:
            case MAP2GRID_DIST_NGP:
            case MAP2GRID_DIST_TSC:
               fwrite(&mass_assignment_scheme,sizeof(int),1,fp_out);
               break;
            default:
               SID_trap_error("Unknown mass assignment scheme (%d) in write_grid().",ERROR_LOGIC,mass_assignment_scheme);
         }
      }
      else
         fp_out=fopen(filename_out,"a");
      fwrite(grid_identifier,sizeof(char),GRID_IDENTIFIER_SIZE,fp_out);
   }

   // Append the grid to the file
   void   *buffer=NULL;
   int     i_rank;
   size_t  alloc_size;
   size_t  alloc_size_local;
   alloc_size_local=(size_t)field->n_field_R_local*sizeof(fftw_real);
   SID_Allreduce(&alloc_size_local,&alloc_size,1,SID_SIZE_T,SID_MAX,SID.COMM_WORLD);
   buffer=SID_malloc(alloc_size);
   SID_set_verbosity(SID_SET_VERBOSITY_RELATIVE,-1);
   //remove_buffer_FFT_R(field);
   SID_set_verbosity(SID_SET_VERBOSITY_DEFAULT);
   for(i_rank=0;i_rank<SID.n_proc;i_rank++){
      size_t alloc_size_i;
      memcpy(buffer,field->field_local,alloc_size_local);
      alloc_size_i=alloc_size_local;
      SID_Bcast(&alloc_size_i,sizeof(size_t),i_rank,SID.COMM_WORLD);
      SID_Bcast(buffer,alloc_size_i,i_rank,SID.COMM_WORLD);
      if(SID.I_am_Master)
         fwrite(buffer,alloc_size_i,1,fp_out);
   }
   SID_free(SID_FARG buffer);
   if(SID.I_am_Master)
      fclose(fp_out);
   SID_log("Done.",SID_LOG_CLOSE);
}

