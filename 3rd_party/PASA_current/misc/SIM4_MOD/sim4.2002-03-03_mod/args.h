#ifndef SIM_ARGS_H
#define SIM_ARGS_H
/* $Id: args.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

typedef struct argv_scores {
	double E;
	int I;
	int M;
	int O;
	int V;
} argv_scores_t;

bool get_argval(int, int *);
bool get_fargval(int, double *);
bool get_strargval(int, char **);
bool get_cargval(int, char **);
void ckargs(const char *, int , char **, int );
void fprintf_argv(FILE* fp);
void ck_argc(const char *);

extern char *argv0;
#endif
