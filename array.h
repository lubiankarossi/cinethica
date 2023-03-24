#ifndef TA_ARRAY
#define TA_ARRAY

#include <stdlib.h>

#define ta_array(type) struct { size_t v_size; type *ptr; }
#define ta_initialize(type,val,n) ((val).v_size = n, (val).ptr = (type*)malloc( sizeof(type) * n ) )
#define ta_get(val, i) ((val).ptr[(i)])
#define ta_destroy(val) free((val).ptr)
#define ta_size(val) ((val).v_size)

#endif // TA_ARRAY
