
#include "libc.h"
#include "types.h"
#include "misc.h"
#include "seq.h"
#include "encoding.h"

static const char rcsid[] =
"$Id: seq.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";

static int parse_fname(const char* arg, 
		char **fname, int *from, int *len, char **maskfile)
{
        char *p = 0;
        int flags = 0;

	/* "seqfile{maskfile}[from,to]-" */
        *fname = copy_string(arg);

        p = (*fname)+strlen(*fname)-1;
        if (*p == '-') {
            *p = 0;
            flags |= SEQ_DO_REVCOMP;
        }

        if ((p = strchr(*fname, '['))) {
                int to;

		if (sscanf(p+1, "%d,%d", from, &to) != 2)
			return -1;
		if (*from <= 0 || *from > to)
			return -1;
		*p = '\0';
		*len = to - *from + 1;
		flags |= (SEQ_DO_SUBRANGE|SEQ_IS_SUBRANGE);
	} else {
		*from = 1;
		*len = -1;
	}

        if ((p = strchr(*fname, '{'))) {
		char *q = strchr(p+1, '}');
		if (q) {
			*p = *q = 0;
			if (maskfile) {
				*maskfile = copy_string(p+1);
				flags |= SEQ_DO_MASK;
			}
		}
	} else {
		*maskfile = copy_string(""); /* XXX ugh */
	}
	return flags;
}

static int check_flags(int flags)
{
	switch (flags & (SEQ_DISALLOW_AMB|SEQ_ALLOW_AMB)) {
	case 0: 
		/* default is to allow ambiguious */
		flags |= SEQ_ALLOW_AMB;
		break;
	case SEQ_ALLOW_AMB:
	case SEQ_DISALLOW_AMB:
		break;
	case SEQ_DISALLOW_AMB|SEQ_ALLOW_AMB:
		fatalf("seq_open: contradictory flags: SEQ_DISALLOW_AMB|SEQ_ALLOW_AMB");
	}
	return flags;
}

SEQ* seq_open(const char *fname, const char *mode, int flags)
{
	SEQ *s = ckallocz(sizeof(SEQ));
	int r;

	mode = 0;

	r = parse_fname(fname, 
		&(s->fname), &(s->from), &(s->slen), &(s->maskname));
	if (r == -1)
		fatalf("improper positions specification: %s", fname);

	s->flags = check_flags(r|flags);
	s->fp = ckopen(s->fname, "r");
	s->count = 0;
	s->offset = 0;
	return s;
}

SEQ *seq_copy(const SEQ *s)
{
	SEQ *ss = ckallocz(sizeof(SEQ));
	*ss = *s;
	ss->seq = (uchar*)copy_string((const char*)s->seq);
	ss->header = copy_string(s->header);
	ss->fname = copy_string(s->fname);
	ss->maskname = copy_string(s->fname);
	ss->fp = 0; /* XXX - no subsequent seq_read operations allowed */
	ss->offset = 0; /* XXX - no subsequent seq_read operations allowed */
	/* alternatively,
		ss->fp = ckopen(ss->fname, "r");
		ckfseek(ss->fp, ckftell(s->fp), SEEK_SET);
	   but this is more expensive.
	*/
	return ss;
}

SEQ *seq_subseq(const SEQ *s, int origin, int length)
{
	SEQ *ss;

	/* XXX - probably should do reference counting. */
	/*       1-indexing is ugly */
	if (origin < 1 || length < 0)
		return 0;
	if (SEQ_LEN(s) < origin+length-1)
		return 0;
	ss = ckallocz(sizeof(SEQ));
	*ss = *s;
	ss->flags = s->flags|SEQ_IS_SUBSEQ;
	ss->fp = 0;
	ss->offset = 0;
	ss->from = 1;
	ss->seq = s->seq + origin - 1; 
	ss->slen = length;
	return ss;
}

SEQ* seq_from_chars(unsigned char *chrs, unsigned int len)
{
	SEQ *s = ckallocz(sizeof(SEQ));
	s->fname = 0;
	s->header = 0;
	s->hlen = 0;
	s->seq = chrs;
	s->from = 1;
	s->slen = len;
	s->maskname = 0;
	s->flags = SEQ_IS_SUBSEQ;
	s->count = 0;
	s->fp = 0;
	s->offset = 0;
	return s;
}

const char *seq_set_header(SEQ *s, const char *h)
{
	if (s && s->header) {
		free(s->header);
		s->header = copy_string(h);
		return s->header;
	}
	return 0;
}

SEQ* seq_close(SEQ *s)
{
	if (s) {
		if (!(s->flags & SEQ_IS_SUBSEQ)) {
			if (s->fp) {
				if (s->flags & SEQ_HAS_PIPE)
					ckpclose(s->fp);
				else
					fclose(s->fp);
			}
			if (s->fname) free(s->fname);
			if (s->header) free(s->header);
			if (s->seq) free(s->seq);
			if (s->maskname) free(s->maskname);
		}
		memset(s, 0, sizeof(SEQ));
		free(s);
	}
	return 0;
}

uchar dna_cmpl(uchar ch)
{
	/* XXX - assumes ascii, returns space on error. */
	return dna_complement[ch];
}

static SEQ *seq_revcomp_helper(SEQ *seq)
{
	uchar *s, *p;

	/* assert(SEQ_CHARS in dcomp-' '); */
	/* seq_read should check this. */

	s = SEQ_CHARS(seq);
	p = s+SEQ_LEN(seq)-1;
	while (s<=p) {
		uchar c;

		c = dna_cmpl(*s); 
		*s = dna_cmpl(*p); 
		*p = c;
		++s, --p;
	}
	return seq;
}

SEQ *seq_revcomp_inplace(SEQ *seq) 
{
	seq_revcomp_helper(seq);
	seq->flags ^= SEQ_IS_REVCOMP;
	return seq;
}

SEQ *seq_get(const char *fname, const char *mode, int flags)
{
	SEQ *s = seq_open(fname, mode, flags);
	int r = seq_read(s);
	if (r < 0)
		fatalfr("could not read from %s", fname);
	else if (r == 0)
		return 0;
	else
		return s;
	/*NOTREACHED*/
	return 0;
}

int seq_count(SEQ *s)
{
	return s->count;
}

int seq_revisit(SEQ *s, long offset)
{
	return fseek(s->fp, offset, SEEK_SET);
}

long seq_offset(SEQ *s)
{
	return s->offset;
}

#ifdef TEST
int main()
{
	printf("%d\n", sizeof(dna_complement));
	exit(0);
}
#endif
