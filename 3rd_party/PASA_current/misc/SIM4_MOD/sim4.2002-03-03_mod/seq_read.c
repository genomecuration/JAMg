#include "libc.h"
#include "types.h"
#include "misc.h"
#include "seq.h"
#include "encoding.h"
#include "charvec.h"

#ifndef __lint
static const char rcsid[] =
"$Id: seq_read.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

static SEQ *seq_mask_inplace(SEQ *seq);
static int getpair(FILE *fp, int *a, int *b);
static char *byte_fill_range(uchar *p, int l, int c, int a, int b);
static int ws(int c);
static int getnwc(FILE *fp);
static void un_getc(int c, FILE *fp);
static void char_append(charvec_t *s, int c);

static int ws(int c)
{
	return (c == ' ' || c == '\t' || c == '\n' || c == '\f' || c == '\r');
}

static int getnwc(FILE *fp)
{
	int c = EOF;
	if (!feof(fp)) do c = getc(fp); while (c != EOF && ws(c));
	return c;
}

static void un_getc(int c, FILE *fp)
{
	if (c != EOF)
		if (ungetc(c, fp) == EOF)
			fatalf("cannot ungetc '%c'", c);
}

static void char_append(charvec_t *s, int c)
{
	if (!charvec_append(s, c))
		fatal("cannot append");
}

int seq_read(SEQ *seq)
{
	int b, c;
	charvec_t shdr = charvec_INIT(ckrealloc,free);
	charvec_t sseq = charvec_INIT(ckrealloc,free);

	if (feof(seq->fp))
		return 0;

	if (seq->count > 0) { 
		if (seq->flags & SEQ_IS_SUBRANGE) {
			return 0;
		} else {
			seq->from = 1;
			seq->slen = -1; /* will be computed below */
		}
	}

	if (seq->header) ZFREE(seq->header);
	if (seq->seq) ZFREE(seq->seq);

	seq->offset = ftell(seq->fp);

	/* --- header --- */
	c = getnwc(seq->fp);
	if (c == '>') {
		while (c != '\n' && c != EOF) {
			char_append(&shdr, c);
			c = getc(seq->fp);
		}
	} else {
		un_getc(c, seq->fp);
	}
	char_append(&shdr, 0);
	seq->header = shdr.a;
	seq->hlen = shdr.len;

	/* --- seq --- */
	b = '\n';
	c = getnwc(seq->fp);
	while ((c != EOF) && !(b == '\n' && c == '>')) {
		switch (nfasta_ctype[c]) {
		case Nfasta_nt:
			char_append(&sseq, c);
			break;
		case Nfasta_ws:
			/* skip space */
			break;
		case Nfasta_amb:
			if (seq->flags & SEQ_ALLOW_AMB) {
				char_append(&sseq, c);
				break;
			}
			/* FALLTHRU */
		default:
			fatalf("non-DNA character '%c' in sequence '%s'", c, seq->fname);
			break;
		}
		b = c;
		c = getc(seq->fp);
	}
	un_getc(c, seq->fp);

	/* check conformance */
	if (SEQ_LEN(seq) == -1) {
		char_append(&sseq, 0);
		charvec_fit(&sseq);
		seq->seq = (uchar*)sseq.a;
		seq->slen = sseq.len;
		if (seq->slen > 0) --seq->slen;  /* don't include '\0' */
	} else {
		charvec_t ssub = charvec_INIT(ckrealloc,free);
		int i;

		if (SEQ_FROM(seq) < 1 ||
		    (int)sseq.len < SEQ_FROM(seq) ||
		    SEQ_TO(seq) < 1 ||
		    (int)sseq.len < SEQ_TO(seq) ||
		    SEQ_TO(seq) < SEQ_FROM(seq))
			fatalf("range [%d,%d] incommensurate with sequence [%d,%d]",
				SEQ_FROM(seq), SEQ_TO(seq), 1, sseq.len);
		
		for (i = SEQ_FROM(seq); i <= SEQ_TO(seq); ++i)
			char_append(&ssub, sseq.a[i-1]);
		char_append(&ssub, 0);
		charvec_fini(&sseq);
		seq->seq = (uchar*)ssub.a;
	}

	seq->flags = seq->flags &~ SEQ_IS_REVCOMP;
        if (seq->flags & SEQ_DO_REVCOMP) {
	    (void)seq_revcomp_inplace(seq);
	}
	if (seq->flags & SEQ_HAS_MASK) {
	    (void)seq_mask_inplace(seq);
	}

	seq->count++;
	return 1;
}


static SEQ *seq_mask_inplace(SEQ *seq)
{
        int a, b;  
        FILE *fp = fopen(seq->maskname, "r");

        if (fp == 0) {
                fatalfr("cannot open '%s'", seq->maskname);
                return 0;
        } else {
                while (getpair(fp, &a, &b))
                        byte_fill_range(SEQ_CHARS(seq),SEQ_LEN(seq),'X',a,b);
                fclose(fp);
                seq->flags |= SEQ_HAS_MASK;
                return seq;
        }
}

#define BUF 128
static int getpair(FILE *fp, int *a, int *b)
{
        char buf[BUF];

        /* XXX - should handle comments, etc */
        if (fgets(buf, (int)sizeof(buf), fp) == 0)
                return 0;
        if (sscanf(buf, "%d%d", a, b) != 2)
                return 0;
        return 1;
}

static char *byte_fill_range(uchar *p, int l, int c, int a, int b)
{
        /* fill [a,b] (1-indexed) in p with c */

        a--; b = b-a; /* make it into a 0-indexed, open interval */
        return (b < 0 || l < b) ? 0 : memset(p+a, c, (size_t)b);
}
