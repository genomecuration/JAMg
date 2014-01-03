#ifndef SIM_SEQ_H
#define SIM_SEQ_H
/* $Id: seq.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

typedef struct seq_data {
	uchar *seq;
	int slen; /* bytes in seq, not including '\0' */
	int origin;
} seq_data_t;

typedef struct seq_file {
	FILE *fp;
	int flags;
	int count;  /* how many contigs we have already read */
	long offset; /* the starting offset of the contig we just read */

	char *maskname;

	char *fname;
	int from; /* 1 based */

	char *header;
	int hlen; /* bytes in header, not including '\0' */

	uchar *seq;
	int slen; /* bytes in seq, not including '\0' */
} SEQ;

#define SEQ_NAME(s) ((s)->fname)
#define SEQ_LEN(s) ((s)->slen)
#define SEQ_TO(s)  ((s)->slen + (s)->from - 1)
#define SEQ_FROM(s) ((s)->from)
#define SEQ_AT(s,i) ((s)->seq[i])
/* #define SEQ_LEN(s) ((s)->to - (s)->from + 1) */

#define SEQ_HEAD(s) ((s)->header)
#define SEQ_HLEN(s) ((s)->hlen)

#define SEQ_CHARS(s) ((s)->seq)
#define SEQ_SAME(a,b) ((a)==(b))

enum /* powerset */
{ SEQ_IS_SUBRANGE = (1<<0) /* seq is a subrange of a file */
, SEQ_IS_REVCOMP  = (1<<1) /* seq is reverse compliment of a file */
, SEQ_IS_SUBSEQ   = (1<<2) /* seq is a reference to another seq */
, SEQ_HAS_MASK    = (1<<3) /* seq has a mask applied */
, SEQ_HAS_PIPE	  = (1<<4) /* input fd is a pipe */
, SEQ_DO_REVCOMP  = (1<<5) /* make it so after open */
, SEQ_DO_SUBRANGE = (1<<6) /* make it so after open */
, SEQ_DO_MASK	  = (1<<7) /* make it so after open */
, SEQ_ALLOW_AMB	  = (1<<8) /* checked while reading */
, SEQ_DISALLOW_AMB  = (1<<9) /* checked while reading */
};

SEQ* seq_open(const char *fname, const char *mode, int flags);
SEQ* seq_close(SEQ *s);
int seq_read(SEQ *seq);
const char *seq_set_header(SEQ *s, const char *h);
SEQ *seq_copy(const SEQ *s);
SEQ *seq_subseq(const SEQ *s, int origin, int length);
SEQ *seq_revcomp_inplace(SEQ *seq);
SEQ *seq_get(const char *fname, const char *mode, int flags);
SEQ* seq_from_chars(unsigned char *chrs, unsigned int len);
uchar dna_cmpl(unsigned char);

int seq_count(SEQ *s);
int seq_revisit(SEQ *s, long offset);
long seq_offset(SEQ *s);

#endif 
