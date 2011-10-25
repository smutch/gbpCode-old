#include <stdlib.h>
#include <string.h>
#include <gbpCommon.h>
#include <gbpSID.h>
#include <gbpDomain.h>

void exchange_ring_buffer(void     *send_buffer,
                          size_t    buffer_type_size,
                          size_t    send_count,
                          void     *receive_buffer,
                          size_t   *receive_count_in,
                          int       i_rank){
  size_t  send_buffer_size;
  size_t  receive_buffer_size;
  size_t  receive_count_temp;
  size_t *receive_count;
  int     rank_to,rank_from;

  // At times, we may not want to bother
  //   returning the number of received items
  //   (perhaps because it's fixed and known)
  //   This takes care of that if receive_count_in
  //   is set to NULL.
  if(receive_count_in==NULL)
    receive_count=&receive_count_temp;
  else
    receive_count=receive_count_in;
#ifdef USE_MPI
  // Exchange buffer sizes
  send_buffer_size=send_count*buffer_type_size;
  if(send_buffer==NULL || send_buffer_size<=0)
    send_buffer_size=0;
  if(i_rank!=0){
    rank_to=SID.My_rank-i_rank;
    if(rank_to<0)                rank_to+=SID.n_proc;
    else if(rank_to>=SID.n_proc) rank_to-=SID.n_proc;
    rank_from=SID.My_rank+i_rank;
    if(rank_from<0)                rank_from+=SID.n_proc;
    else if(rank_from>=SID.n_proc) rank_from-=SID.n_proc;
    MPI_Sendrecv(&send_count,
                 1,
                 MPI_SIZE_T,
                 rank_to,
                 1236269,
                 receive_count,
                 1,
                 MPI_SIZE_T,
                 rank_from,
                 1236269,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    receive_buffer_size=(*receive_count)*buffer_type_size;
    MPI_Sendrecv(send_buffer,
                 send_buffer_size,
                 MPI_BYTE,
                 rank_to,
                 1256269,
                 receive_buffer, 
                 receive_buffer_size,
                 MPI_BYTE,
                 rank_from,
                 1256269,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
  }
  else{
    receive_buffer_size=send_buffer_size;
    (*receive_count)   =send_count;
    memcpy(receive_buffer,send_buffer,receive_buffer_size);
  }
#else
  if(send_buffer!=NULL && send_buffer_size>0){
    receive_buffer_size=send_buffer_size;
    (*receive_count)   =send_count;
    memcpy(receive_buffer,send_buffer,receive_buffer_size);
  }
#endif
}

