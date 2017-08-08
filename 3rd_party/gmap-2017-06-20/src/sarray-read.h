/* $Id: sarray-read.h 207324 2017-06-14 19:41:18Z twu $ */
#ifndef SARRAY_READ_INCLUDED
#define SARRAY_READ_INCLUDED
#include "access.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"


#define T Sarray_T
typedef struct T *T;

/* For benchmarking */
Univcoord_T
Sarray_size (Sarray_T this);

#if 0
extern void
Sarray_shmem_remove (char *dir, char *fileroot, char *snps_root, Mode_T mode, bool fwdp);
#endif

extern Univcoord_T *
Sarray_array (T this);

extern Univcoord_T
Sarray_position (T sarray, Sarrayptr_T i);

extern T
Sarray_new (char *dir, char *fileroot, Access_mode_T sarray_access, Access_mode_T lcp_access,
	    Access_mode_T guideexc_access, Access_mode_T indexij_access, bool sharedp, Mode_T mode, bool fwdp);
extern void
Sarray_free (T *old);

extern void
Sarray_read (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, bool *successp,
	     UINT4 *nmatches, char *query, UINT4 querylength, int queryoffset,
	     Compress_T query_compress, T sarray, bool plusp, int genestrand,
	     char conversion[]);

extern Univcoord_T *
Sarray_lookup (int *nhits, T sarray, char *query, UINT4 querylength, int queryoffset,
	       Compress_T query_compress, bool plusp, int genestrand,
	       char conversion[]);

#undef T
#endif


