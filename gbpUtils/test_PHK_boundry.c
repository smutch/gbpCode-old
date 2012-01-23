#define  _MAIN
#include <gbpLib.h>
#include <gbpPHKs.h>

int main(int argc,char *argv[]){
   int n_bits;
   int key_start;
   int key_stop;

   SID_init(&argc,&argv,NULL);
 
   // Read inputs
   n_bits   =(int)atoi(argv[1]);
   key_start=(int)atoi(argv[2]);
   key_stop =(int)atoi(argv[3]);
   SID_log("Writing %d-bit Peano-Hilbert key information for range %d->%d...",SID_LOG_OPEN,n_bits,key_start,key_stop);
   SID_log("Maximum dimensional size: %d",SID_LOG_COMMENT,PHK_DIM_SIZE(n_bits));
   SID_log("Maximum key value:        %d",SID_LOG_COMMENT,PHK_DIM_SIZE(n_bits)*PHK_DIM_SIZE(n_bits)*PHK_DIM_SIZE(n_bits));

   // Write the positions of the keys
   PHK_t  key;
   FILE  *fp;
   int    i_x;
   int    i_y;
   int    i_z;
   
   SID_log("Writing domain keys to file {keys_domain.dat}...",SID_LOG_OPEN);
   fp=fopen("keys_domain.dat","w");
   fprintf(fp,"# Peano-Hilbert domain key information\n");
   fprintf(fp,"#   No. of bits/dim =%d\n",    n_bits);
   fprintf(fp,"#   Range           =%d->%d\n",key_start,key_stop);
   fprintf(fp,"#   Dimensional size=%d\n",    PHK_DIM_SIZE(n_bits));
   fprintf(fp,"#   Max. key value  =%d\n",    PHK_DIM_SIZE(n_bits)*PHK_DIM_SIZE(n_bits)*PHK_DIM_SIZE(n_bits));
   fprintf(fp,"#\n");
   fprintf(fp,"#   Columns: (1) key\n");
   fprintf(fp,"#            (2) x-position\n");
   fprintf(fp,"#            (3) y-position\n");
   fprintf(fp,"#            (4) z-position\n");
   for(key=key_start;key<=key_stop;key++){
      compute_PHK_to_Cartesian(n_bits,key,&i_x,&i_y,&i_z);
      fprintf(fp,"%7zd %7d %7d %7d\n",key,i_x,i_y,i_z);
   }
   fclose(fp);
   SID_log("Done.",SID_LOG_CLOSE);

   // Write the positions of the boundary keys
   PHK_t *keys_boundary;
   int    n_keys_boundary;
   int    i_key;
   SID_log("Writing boundary keys to file {keys_boundary.dat}...",SID_LOG_OPEN);
   compute_PHK_boundary_keys(n_bits,key_start,key_stop,&n_keys_boundary,&keys_boundary);
   fp=fopen("keys_boundary.dat","w");
   fprintf(fp,"# Peano-Hilbert boundary key information\n");
   fprintf(fp,"#   No. of bits/dim =%d\n",    n_bits);
   fprintf(fp,"#   Range           =%d->%d\n",key_start,key_stop);
   fprintf(fp,"#   Dimensional size=%d\n",    PHK_DIM_SIZE(n_bits));
   fprintf(fp,"#   Max. key value  =%d\n",    PHK_DIM_SIZE(n_bits)*PHK_DIM_SIZE(n_bits)*PHK_DIM_SIZE(n_bits));
   fprintf(fp,"#\n");
   fprintf(fp,"#   Columns: (1) key\n");
   fprintf(fp,"#            (2) x-position\n");
   fprintf(fp,"#            (3) y-position\n");
   fprintf(fp,"#            (4) z-position\n");
   for(i_key=0;i_key<n_keys_boundary;i_key++){
      compute_PHK_to_Cartesian(n_bits,keys_boundary[i_key],&i_x,&i_y,&i_z);
      fprintf(fp,"%7zd %7d %7d %7d\n",keys_boundary[i_key],i_x,i_y,i_z);
   }
   fclose(fp);
   SID_free(SID_FARG keys_boundary);
   SID_log("Done.",SID_LOG_CLOSE);
   SID_log("Done.",SID_LOG_CLOSE);

   SID_exit(ERROR_NONE);
}

