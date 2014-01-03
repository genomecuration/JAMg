#include "libc.h"
#include "types.h"
#include "args.h"
#include "seq.h"
#include "dna.h"
#include "misc.h"
#include "prnt.h"

#ifndef __lint
static const char rcsid[] =
"$Id: prnt.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

/* XXX */
static int offset1;
static int offset2;

enum { BUFSIZE=128 };
#define ckprintf (void)printf /* XXX */

static char *subseq_label(char *buf, unsigned int size, int n);
static const char* revflag(SEQ *s);
static const char* revlabel(SEQ *s);

static void print_align_header_n(SEQ *seq1, SEQ *seq2, argv_scores_t *ds, int n)
{
        int f, t, F, T;
        char buf[BUFSIZE];

        ckprintf("#:lav\n\nd {\n  \"");
        ck_argc("print_align_header");
        fprintf_argv(stdout);
        ckprintf("\n   M = %d, I = %d, V = %d", ds->M, ds->I, ds->V);
        ckprintf(", O = %d, E = %g", ds->O, ds->E);
        ckprintf("\"\n}\n");

        if (get_argval('f', &f)) {
                if (!get_argval('t',&t) || !get_argval('F',&F) ||
                    !get_argval('T',&T))
                        fatal("Inconsistent use of `f`, `t`, `F', `T' args.");
                offset1 = SEQ_FROM(seq1) - f;
                offset2 = SEQ_FROM(seq2) - F;
        } else {
                f = SEQ_FROM(seq1); t = SEQ_TO(seq1);
                F = SEQ_FROM(seq2); T = SEQ_TO(seq2);
                offset1 = offset2 = 0;
        }
        ckprintf("s {\n  \"%s%s\" %d %d\n  \"%s%s\" %d %d\n}\n",
                        SEQ_NAME(seq1), revflag(seq1), f, t,
                        SEQ_NAME(seq2), revflag(seq2), F, T);
        ckprintf("h {\n   \"%s%s\"\n   \"%s%s%s\"\n}\n",
                SEQ_HEAD(seq1),
                revlabel(seq1),
                SEQ_HEAD(seq2),
                revlabel(seq2),
                subseq_label(buf, sizeof buf, n));
}

/* print_align_header  -------------  print the top part of an alignment file */
void print_align_header(SEQ *seq1, SEQ *seq2, argv_scores_t *ds)
{
	print_align_header_n(seq1, seq2, ds, 0);
}

static char *subseq_label(char *buf, unsigned int size, int n)
{
	assert(size > 0);
	buf[0] = 0;
	if (n > 0) snprintf(buf, size, " (subsequence #%d)", n);
	return buf;
}

static const char* revflag(SEQ *s)
{
	return (s->flags & SEQ_IS_REVCOMP) ? "-" : "";
}

static const char* revlabel(SEQ *s)
{
	return (s->flags & SEQ_IS_REVCOMP) ? " (reverse complement)" : "";
}


/* print_align  ----------------------------------- print a general alignment */
void print_align(int score, uchar *seq1, uchar *seq2, int beg1, int end1, int beg2, int end2,int *S)
{
        int M, N, i, j, op, start_i, start_j, match, run, pct;
        uchar *P, *p, *q;

        beg1 += offset1;
        end1 += offset1;
        beg2 += offset2;
        end2 += offset2;

        M = end1 - beg1 + 1;
        N = end2 - beg2 + 1;
        ckprintf("a {\n  s %d\n  b %d %d\n  e %d %d\n",
                score, beg1, beg2, end1, end2);
        for (i = j = 0; i < M || j < N; ) {
                start_i = i;
                start_j = j;
                match = 0;
                P= p = seq1 + beg1 + i - 1;
                q = seq2 + beg2 + j - 1;
                while (i < M && j < N && *S == 0) {
                        if (*p++ == *q++)
                                ++match;
                        ++i;
                        ++j;
                        ++S;
                }
                run = p - P;
                pct = (run > 0) ? ((100*match + run/2)/run) : 0; /* round */
                ckprintf("  l %d %d %d %d %d\n",
                        beg1+start_i, beg2+start_j, beg1+i-1, beg2+j-1, pct);
                if (i < M || j < N) {
                        if ((op = *S++) > 0) j += op; else i -= op;
                }       
        }       
        ckprintf("}\n");
}
