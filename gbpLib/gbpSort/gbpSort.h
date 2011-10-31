#define SORT_LOCAL               1
#define SORT_GLOBAL              0
#define SORT_COMPUTE_INPLACE     1
#define SORT_COMPUTE_NOT_INPLACE 0
#define SORT_COMPUTE_RANK        100
#define SORT_COMPUTE_INDEX       101
#define SORT_INPLACE_ONLY        102

void sort(void          *sval,
	  size_t         nval,
	  size_t       **index,
	  SID_Datatype   type,
	  int            flag_local,
	  int            flag_compute_index,
	  int            flag_in_place);
void heap_sort(void    *data_in,
               size_t   n_data,
               size_t **index,
               SID_Datatype      data_type,
               int      flag_compute_index,
               int      flag_in_place);
void merge_helper(void   *data,
                  size_t *index,
                  size_t  left,
                  size_t  right,
                  void   *scratch_d,
                  size_t *scratch_i,
                  SID_Datatype     data_type,
                  int     flag_compute_index);
void merge_sort(void    *data_in,
                size_t   n_data,
                size_t **index,
                SID_Datatype      data_type,
                int      flag_compute_index,
                int      flag_in_place);
