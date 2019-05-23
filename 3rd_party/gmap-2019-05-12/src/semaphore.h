/* $Id: semaphore.h 183951 2016-02-09 01:00:57Z twu $ */
#ifndef SEMAPHORE_INCLUDED
#define SEMAPHORE_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>		/* For key_t and ushort */


#define NSEMAPHORES 2
#define SEMNO_NA 0
#define SEMNO_LOCK 0
#define SEMNO_KEEP 1

#define SEMAPHORE_RESIDENT 99
#define SEMAPHORE_FREEABLE 1


extern int
Semaphore_create (key_t key, int nsems, ushort *vals);

extern int
Semaphore_find (key_t key);

extern void
Semaphore_delete (int semid);

extern void
Semaphore_set_value (int semid, int semno, int value);

extern int
Semaphore_get_value (int semid, int semno);

extern int
Semaphore_nwaiting (int semid, int semno);

extern void
Semaphore_lock (int semid, int semno);

extern void
Semaphore_unlock (int semid, int semno);

#endif

