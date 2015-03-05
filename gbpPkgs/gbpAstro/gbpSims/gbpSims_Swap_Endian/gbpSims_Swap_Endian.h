#ifndef GBPSSIMPL_AWAKE
#define GBPSSIMPL_AWAKE

#define SWAP_SSIMPL_ENDIAN_TO_NATIVE   TTTP00
#define SWAP_SSIMPL_ENDIAN_FROM_NATIVE TTTP01

// Function definitions
#ifdef __cplusplus
extern "C" {
#endif
int  swap_endian_snapshot(const char *filename_in,const char *filename_out,int region,int snap_number,int mode,int *IDs_byte_size);
int  swap_endian_smooth(const char *filename_in,const char *filename_out,int region,int snap_number,int mode,int IDs_byte_size);
void swap_endian_grids(const char *filename_in,const char *filename_out,int mode);
void swap_endian_halos(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode);
void swap_endian_catalogs(const char *filename_in_root,const char *filename_out_root,const char *filename_halo_type,int snap_number,int mode);
#ifdef __cplusplus
}
#endif

#endif

