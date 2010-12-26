#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpDomain.h>

void exchange_slab_buffer_right(void      *send_buffer,
                                size_t     send_buffer_size,
                                void      *receive_buffer,
                                size_t    *receive_buffer_size,
                                slab_info *slab){
  size_t size_temp;
#ifdef USE_MPI
  // Exchange buffer sizes
  if(send_buffer==NULL || send_buffer_size<=0)
    size_temp=0;
  else
    size_temp=send_buffer_size;
  MPI_Sendrecv(&size_temp,
               1,
               MPI_SIZE_T,
               slab->rank_to_right,
               123,
               receive_buffer_size,
               1,
               MPI_SIZE_T,
               slab->rank_to_left,
               123,
               MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
  if(size_temp>0 || *receive_buffer_size>0)
    MPI_Sendrecv(send_buffer,
                 send_buffer_size,
                 MPI_BYTE,
                 slab->rank_to_right,
                 125,
                 receive_buffer, 
                 (*receive_buffer_size),
                 MPI_BYTE,
                 slab->rank_to_left,
                 125,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
#else
  if(send_buffer!=NULL && send_buffer_size>0){
    (*receive_buffer_size)=send_buffer_size;
    memcpy(receive_buffer,send_buffer,(*receive_buffer_size));
  }
#endif
}

