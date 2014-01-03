#define _XOPEN_SOURCE /* tell sun we want popen, etc */
#include "libc.h"
#include "types.h"
#include "misc.h"
#include "args.h"

#ifndef __lint
static const char rcsid[] =
"$Id: misc.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

/* fatal ---------------------------------------------- print message and die */
void fatal(const char *msg)
{
	fflush(stdout);
	fatalf("%s", msg);
	exit(1);
}

/* fatalf --------------------------------- format message, print it, and die */

static void print_argv0(void)
{
	if (argv0) {
		char *p = strrchr(argv0, '/');
		(void)fprintf(stderr, "%s: ", p ? p+1 : argv0);
	}
}

void fatalf(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fflush(stdout);
	print_argv0();
	(void)vfprintf(stderr, fmt, ap);
	(void)fputc('\n', stderr);
	va_end(ap);
	exit(1);
}

void fatalfr(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fflush(stdout);
	print_argv0();
	(void)vfprintf(stderr, fmt, ap);
	(void)fprintf(stderr, ": %s\n", strerror(errno));
	va_end(ap);
	exit(1);
}

int psublast_debug = 0;

void debugf(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	if (psublast_debug) {
		fflush(stdout);
		print_argv0();
		if (vfprintf(stderr, fmt, ap) < 0)
			exit(1);
	}
	va_end(ap);
}

void debugff(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	if (psublast_debug) {
		fflush(stdout);
		print_argv0();
		if (vfprintf(stderr, fmt, ap) < 0)
			exit(1);
		if (fflush(stderr) != 0)
			exit(1);
	}
	va_end(ap);
}

/* ckopen -------------------------------------- open file; check for success */
FILE *ckopen(const char *name, const char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)
		fatalfr("Cannot open %s.", name);
	return fp;
}

/* ckalloc -------------------------------- allocate space; check for success */
void *ckalloc(size_t amount)
{
	void *p;

	assert((long)amount >= 0); 
	if (amount == 0)
		amount = 1; /* ANSI portability hack */
	if ((p = malloc(amount)) == NULL)
		fatalf("Ran out of memory trying to allocate %lu.",
			(unsigned long)amount);
#if 0
	memset(p, 0, amount); /* XXX */
#endif
	return p;
}

/* ckallocz -------------------- allocate space; zero fill; check for success */
void *ckallocz(size_t amount)
{
	void *p = ckalloc(amount);
	memset(p, 0, amount);
	return p;
}

void *ckfree(void *p)
{
	free(p);
	return 0;
}

/* strsame --------------------------- tell whether two strings are identical */
bool same_string(const char *s, const char *t)
{
	return (strcmp(s, t) == 0);
}

/* strsave -------------------------- save string s somewhere; return address */
char *copy_string(const char *s)
{
	char *p;

	p = ckalloc(strlen(s)+1);	/* +1 to hold '\0' */
	return strcpy(p, s);
}

char *copy_substring(const char *s, int n)
{
	char *p = ckalloc((size_t)n+1);	/* +1 to hold '\0' */
	memcpy(p, s, (size_t)n);
	p[n] = 0;
	return p;
}

long ckftell(FILE *f)
{
	long r = ftell(f);
	if (r < 0)
		fatalfr("bad ftell: %s");
	return r;
}

int ckfseek(FILE *f, long i, int m)
{
	int r = fseek(f, i, m);
	if (r < 0)
		fatalfr("bad fseek: %s");
	return r;
}

void *ckrealloc(void * p, size_t size)
{
	p = p ? realloc(p, size) : malloc(size);
	if (!p)
		fatal("ckrealloc failed");
	return p;
}

FILE *ckpopen(const char *name, const char *mode)
{
	FILE *fp;

	if ((fp = popen(name, mode)) == NULL)
		fatalfr("Cannot open %s.", name);
	return fp;
}

void ckpclose(FILE *fp)
{
	int r = pclose(fp);
	if (r != 0)
		fatalfr("pclose failed (status %d)", r);
}
