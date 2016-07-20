#ifndef GBPPARSE_PARAMTER_FILES_AWAKE
#define GBPPARSE_PARAMTER_FILES_AWAKE

#include <gbpSID.h>
#include <gbpParse_core.h>

#define PARAMETER_STRING_LENGTH  64
#define PARAMETER_MODE_OPTIONAL  TTTP01
#define PARAMETER_MODE_MANDITORY TTTP02
#define PARAMETER_MODE_UNIQUE    TTTP03
#define PARAMETER_MODE_DEFAULT   PARAMETER_MODE_MANDITORY|PARAMETER_MODE_UNIQUE

typedef struct parameter_item_info parameter_item_info;
struct parameter_item_info {
   char                 name[PARAMETER_STRING_LENGTH];
   SID_Datatype         data_type;
   int                  data_size;
   int                  n_read;
   int                  mode;
   int                  flag_set;
   void                *data;
   parameter_item_info *next;
};

typedef struct parameter_list_info parameter_list_info;
struct parameter_list_info {
   parameter_item_info *first;
   parameter_item_info *last;
   int                  count;
};

// Function declarations
#ifdef __cplusplus
extern "C" {
#endif
int  remove_parameter_character(char *line);
void init_parameter_list(parameter_list_info **param_list);
void init_parameter_item(parameter_item_info **param_item,
                         const char           *name,
                         SID_Datatype          data_type,
                         int                   mode);
void free_parameter_item(parameter_item_info **param_item);
int  free_parameter_list(parameter_list_info **param_list);
void add_parameter_to_list(parameter_list_info *param_list,
                           const char          *name,
                           SID_Datatype         data_type,
                           int                  mode);
int  fetch_parameter_data(parameter_list_info *param_list,const char *name,void *value);
int  fetch_parameter_list_item(parameter_list_info *param_list,const char *name,parameter_item_info **item);
void read_gbpParam_file(const char          *filename_in,
                        parameter_list_info *param_list);
#ifdef __cplusplus
}
#endif

#endif

