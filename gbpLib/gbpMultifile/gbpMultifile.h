#ifndef MULTIFILE_AWAKE
#define MULTIFILE_AWAKE
#include <gbpCommon.h>
#include <gbpSID.h>

// V Preprocessor definitions V
// A Preprocessor definitions A

// V --- Datatype definitions --- V
// This datastructure describes the multifile file-pointer
typedef struct fp_multifile_info fp_multifile_info;
struct fp_multifile_info{
   char    filename_root[MAX_FILENAME_LENGTH];
   char    filename_base[MAX_FILENAME_LENGTH];
   FILE   *fp_multifile;
   size_t  data_size;
   int     i_file;
   int     n_files;
   int     n_items_total;
   int     i_item;
   int     i_item_start;
   int     i_item_stop;
   int     n_items_file;
   int     flag_multifile;
};
// A --- Datatype definitions --- A

#ifdef __cplusplus
extern "C" {
#endif   
// V --- ANSI-C function definitions --- V
int fopen_multifile(const char        *filename_multifile_root,
                    size_t             data_size,
                    fp_multifile_info *fp_out,
                    ...);
int fopen_multifile_nth_file(fp_multifile_info *fp_in,
                             int                n);
void fclose_multifile(fp_multifile_info *fp_in);
int fread_multifile(fp_multifile_info *fp_in,
                    void              *data_out,
                    int                item_index);
// A --- ANSI-C function definitions --- A
#ifdef __cplusplus
}
#endif

#endif
