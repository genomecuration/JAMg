static char rcsid[] = "$Id: semaphore.c 183989 2016-02-09 18:07:52Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "semaphore.h"
#include <stdlib.h>
#include <stdio.h>

#include <sys/ipc.h>
#include <sys/sem.h>		/* For semaphores */


int
Semaphore_create (key_t key, int nsems, ushort *values) {
  int semid;
  union semun {
    int val;
    struct semid_ds *buf;
    ushort *array;
  } arg;

  if ((semid = semget(key,nsems,IPC_CREAT | IPC_EXCL | 0666)) == -1) {
    return -1;
  } else {
    arg.array = values;
    semctl(semid,/*semno*/0,SETALL,arg);
    return semid;
  }
}

int
Semaphore_find (key_t key) {
  return semget(key,/*nsems*/0,/*flags*/0);
}


void
Semaphore_delete (int semid) {
  if (semctl(semid,0,IPC_RMID,NULL) == -1) {
    fprintf(stderr,"Error releasing semaphore %d\n",semid);
    /* abort(); */
  }
  return;
}


void
Semaphore_set_value (int semid, int semno, int value) {
  union semun {
    int val;
    struct semid_ds *buf;
    ushort *array;
  } arg;

  arg.val = value;
  semctl(semid,semno,SETVAL,arg);
  return;
}

int
Semaphore_get_value (int semid, int semno) {
  /* printf("Value of semno %d is %d\n",semno,semctl(semid,semno,GETVAL,NULL)); */
  return semctl(semid,semno,GETVAL,/*arg*/NULL);
}


int
Semaphore_nwaiting (int semid, int semno) {
  printf("nwaiting = %d\n",semctl(semid,semno,GETNCNT,NULL));
  return semctl(semid,semno,GETNCNT,/*arg*/NULL);
}


/* If already locked, then puts process to sleep */
void
Semaphore_lock (int semid, int semno) {
  struct sembuf op;

  /* printf("Process %d locking semaphore %d\n",getpid(),semid); */
  op.sem_num = semno;
  op.sem_op = -1;
  op.sem_flg = SEM_UNDO;
  semop(semid,&op,1);

  return;
}

void
Semaphore_unlock (int semid, int semno) {
  struct sembuf op;

  /* printf("Process %d unlocking semaphore %d\n",getpid(),semid); */
  op.sem_num = semno;
  op.sem_op = +1;
  op.sem_flg = SEM_UNDO;
  semop(semid,&op,1);

  /* return semctl(semid,SEMAPHORE_CREATION,GETNCNT,NULL); */
  return;
  /* printf("%d processes still waiting\n",semctl(semid,SEMAPHORE_CREATION,GETNCNT,NULL)); */
}


