#include <gbpCommon.h>
#include <gbpSID.h>

void SID_Type_size(SID_Datatype type,int *size){
  #if USE_MPI
  MPI_Type_size(type,size);
  #else
  if(type==SID_INT)
    (*size)=sizeof(int);
  else if(type==SID_LONG_LONG)
    (*size)=sizeof(long long);
  else if(type==SID_UNSIGNED)
    (*size)=sizeof(unsigned int);
  else if(type==SID_SIZE_T)
    (*size)=sizeof(size_t);
  else if(type==SID_FLOAT)
    (*size)=sizeof(float);
  else if(type==SID_DOUBLE)
    (*size)=sizeof(double);
  else if(type==SID_CHAR)
    (*size)=sizeof(char);
  else
    SID_trap_error("Unsupported SID_Datatype (%d) in SID_Type_size().",ERROR_LOGIC,type);
  #endif
}
