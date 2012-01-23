// Free memory allocated by initialize_field()
#include <gbpLib.h>
#include <gbpPHKs.h>

void compute_PHK_to_Cartesian(int PHK_bit_size,PHK_t key,int *i_x,int *i_y,int *i_z){
   hikey_t  basecoord[3];

   hilbert_i2c(3,(hikey_t)PHK_bit_size,(hikey_t)key,basecoord);
   (*i_x)=(int)basecoord[0];
   (*i_y)=(int)basecoord[1];
   (*i_z)=(int)basecoord[2];

}

