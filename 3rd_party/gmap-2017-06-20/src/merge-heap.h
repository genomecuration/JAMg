#ifndef MERGE_HEAP_INCLUDED
#define MERGE_HEAP_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "list.h"
#include "intlist.h"
#include "merge.h"		/* For Record_T */


extern Record_T *
Merge_records_heap (int *nelts1, List_T stream_list, Intlist_T streamsize_list,
		    Intlist_T querypos_list, Intlist_T diagterm_list, 
		    struct Record_T *all_records);

#endif

