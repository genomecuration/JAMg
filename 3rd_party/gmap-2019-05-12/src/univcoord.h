/* $Id: univcoord.h 218286 2019-01-23 16:46:55Z twu $ */
#ifndef UNIVCOORD_INCLUDED
#define UNIVCOORD_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For SIZEOF_UNSIGNED_LONG_LONG, SIZEOF_UNSIGNED_LONG needed for HAVE_64_BIT */
#endif

#include "types.h"		/* For HAVE_64_BIT */

/* For intervals and IIT files */
#if defined(LARGE_GENOMES) || (defined(UTILITYP) && defined(HAVE_64_BIT))
#include "uint8list.h"
typedef UINT8 Univcoord_T;
typedef Uint8list_T Univcoordlist_T;
#define Univcoordlist_push Uint8list_push
#define Univcoordlist_pop Uint8list_pop
#define Univcoordlist_length Uint8list_length
#define Univcoordlist_reverse Uint8list_reverse
#define Univcoordlist_append Uint8list_append
#define Univcoordlist_next Uint8list_next
#define Univcoordlist_head Uint8list_head
#define Univcoordlist_free Uint8list_free
#define Univcoordlist_keep_one Uint8list_keep_one
#define Univcoordlist_min Uint8list_min
#define Univcoordlist_max Uint8list_max
#define Univcoordlist_to_array Uint8list_to_array
#define Univcoordlist_to_array_out Uint8list_to_array_out
#define Univcoordlist_fill_array Uint8list_fill_array
#define Univcoordlist_fill_array_and_free Uint8list_fill_array_and_free
#define Univcoordlist_to_string Uint8list_to_string
#define Univcoordlist_to_string_offset Uint8list_to_string_offset

#ifdef GSNAP
#include "uint8listpool.h"
typedef Uint8listpool_T Univcoordlistpool_T;
#define Univcoordlistpool_new Uint8listpool_new
#define Univcoordlistpool_free Uint8listpool_free
#define Univcoordlistpool_free_memory Uint8listpool_free_memory
#define Univcoordlistpool_reset Uint8listpool_reset
#define Univcoordlistpool_push Uint8listpool_push
#define Univcoordlistpool_pop Uint8listpool_pop
#define Univcoordlistpool_copy Uint8listpool_copy
#define Univcoordlistpool_copy_but_last Uint8listpool_copy_but_last
#endif

#include "uint8table_rh.h"
typedef Uint8table_T Univcoordtable_T;
#define Univcoordtable_new Uint8table_new
#define Univcoordtable_free Uint8table_free
#define Univcoordtable_put_and_save Uint8table_put_and_save
#define Univcoordtable_get Uint8table_get
#define Univcoordtable_saved_values Uint8table_saved_values
#define Univcoordtable_length Uint8table_length
#define Univcoordtable_keys Uint8table_keys

#else
#include "uintlist.h"
typedef UINT4 Univcoord_T;
typedef Uintlist_T Univcoordlist_T;
#define Univcoordlist_push Uintlist_push
#define Univcoordlist_pop Uintlist_pop
#define Univcoordlist_length Uintlist_length
#define Univcoordlist_reverse Uintlist_reverse
#define Univcoordlist_append Uintlist_append
#define Univcoordlist_next Uintlist_next
#define Univcoordlist_head Uintlist_head
#define Univcoordlist_free Uintlist_free
#define Univcoordlist_keep_one Uintlist_keep_one
#define Univcoordlist_min Uintlist_min
#define Univcoordlist_max Uintlist_max
#define Univcoordlist_to_array Uintlist_to_array
#define Univcoordlist_to_array_out Uintlist_to_array_out
#define Univcoordlist_fill_array Uintlist_fill_array
#define Univcoordlist_fill_array_and_free Uintlist_fill_array_and_free
#define Univcoordlist_to_string Uintlist_to_string
#define Univcoordlist_to_string_offset Uintlist_to_string_offset

#ifdef GSNAP
#include "uintlistpool.h"
typedef Uintlistpool_T Univcoordlistpool_T;
#define Univcoordlistpool_new Uintlistpool_new
#define Univcoordlistpool_free Uintlistpool_free
#define Univcoordlistpool_free_memory Uintlistpool_free_memory
#define Univcoordlistpool_reset Uintlistpool_reset
#define Univcoordlistpool_push Uintlistpool_push
#define Univcoordlistpool_pop Uintlistpool_pop
#define Univcoordlistpool_copy Uintlistpool_copy
#define Univcoordlistpool_copy_but_last Uintlistpool_copy_but_last
#endif

#include "uinttable_rh.h"
typedef Uinttable_T Univcoordtable_T;
#define Univcoordtable_new Uinttable_new
#define Univcoordtable_free Uinttable_free
#define Univcoordtable_put_and_save Uinttable_put_and_save
#define Univcoordtable_get Uinttable_get
#define Univcoordtable_saved_values Uinttable_saved_values
#define Univcoordtable_length Uint8table_length
#define Univcoordtable_keys Uint8table_keys

#endif

#endif

