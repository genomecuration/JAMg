#ifndef SIM_MISC_H
#define SIM_MISC_H
/* $Id: misc.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

#ifdef __GNUC__
#define NORETURN __attribute__((__noreturn__))
#else
#define NORETURN /* */
#endif

#define CLEN(s) (sizeof((s))-1)

/*@exits@*/ void fatal(const char *msg) NORETURN;
/*@exits@*/ void fatalf(const char *fmt, ...) NORETURN;
/*@exits@*/ void fatalfr(const char *fmt, ...) NORETURN;
void debugf(const char *fmt, ...);
void debugff(const char *fmt, ...);

FILE *ckpopen(const char *name, const char *mode);
void ckpclose(FILE*);

FILE *ckopen(const char *name, const char *mode);
/*@only@*/ void *ckalloc(size_t amount);
/*@only@*/ void *ckallocz(size_t amount);
void *ckfree(void *p);
bool same_string(const char *s, const char *t);
char *copy_string(const char *s);
char *copy_substring(const char *s, int n);
long ckftell(FILE*);
int ckfseek(FILE*,long,int);
void *ckrealloc(void *, size_t);

#define ZFREE(p) /*CONSTCOND*/do{free(p);(p)=0;}while(0)

#ifndef RCSID
#define RCSID(id) static const char rcsid[] = id
#endif

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#undef ICEIL
#define ICEIL(x,y) ((((x)-1)/(y))+1)

extern int psublast_debug;

#endif
