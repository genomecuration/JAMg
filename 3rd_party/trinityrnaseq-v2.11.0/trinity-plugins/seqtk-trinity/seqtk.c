/* The MIT License

   Copyright (c) 2008-2016 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)
KHASH_SET_INIT_INT64(64)

typedef kh_reg_t reghash_t;

unsigned read_type = 0; // for trinity in setting /1 or /2


reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	str = calloc(1, sizeof(kstring_t));
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int beg = -1, end = -1;
		reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s);
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s);
						if (end < 0) end = -1;
					}
				}
			}
		}
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg < 0) beg = 0, end = INT_MAX;
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return h;
}

void stk_reg_destroy(reghash_t *h)
{
	khint_t k;
	if (h == 0) return;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

/* constant table */

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

unsigned char seq_nt6_table[256] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

char *seq_nt16_rev_table = "XACMGRSVTWYHKDBN";
unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
unsigned char seq_nt16comp_table[] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static void stk_printstr(const kstring_t *s, unsigned line_len)
{
	if (line_len != UINT_MAX && line_len != 0) {
		int i, rest = s->l;
		for (i = 0; i < s->l; i += line_len, rest -= line_len) {
			putchar('\n');
			if (rest > line_len) fwrite(s->s + i, 1, line_len, stdout);
			else fwrite(s->s + i, 1, rest, stdout);
		}
		putchar('\n');
	} else {
		putchar('\n');
		puts(s->s);
	}
}


static unsigned ctoi(char c) {
    char c2[2] = {'0', '\0'};
    
    c2[0] = c;
    
    unsigned i = atoi(c2);
    
    return(i);
}

static void update_read_name_for_Trinity(char* name, int comment_length, const char* comment) {


    // see if already in the old format:
    //  ie.
    //  61DFRAAXX100204:1:101:11674:10443/1


    int name_len = strlen(name);
    char read_type_char = name[name_len - 1];
    
    // check for _forward or _reverse in read name.

    char* found;
    if ( (found = strstr(name, "_forward")) != NULL) {
        name_len = found - name;
    } else if ( (found = strstr(name, "_reverse")) != NULL) {
        name_len = found - name;
    }
    
    
    if ( found == NULL 
         && 
         name[name_len -2 ] == '/'
         &&
        (read_type_char == '1' || read_type_char == '2') ) {

        // ensure as expected.
        if (ctoi(read_type_char) != read_type) {
            fprintf(stderr, "Error, found read_type %c but expecting read_type %i\n", read_type_char, read_type);
            exit(2);
        } 
        
        return;
    }
    

    // see if in the new format
    //  ie. 
    //  M01581:927:000000000-ARTAL:1:1101:19874:2078 1:N:0:1
    
    read_type_char = comment[0];

    if (found == NULL 
        &&
        comment_length > 1
         &&
             comment[1] == ':'
             &&
             (read_type_char == '1' || read_type_char == '2') ) {
        
        // ensure as expected
        if (ctoi(read_type_char) != read_type) {
            fprintf(stderr, "Error, found read_type %c but expecting read_type %i\n", read_type_char, read_type);
            exit(2);
        } 
        
        // recognized as new format.  Convert to old format that trinity likes.
        name[name_len] = '/';
        name[name_len+1] = comment[0];
        name[name_len + 2] = '\0';
        
        return;
    }

    
        
    // used to error here as per below.
    
    // if got here, none of the above read name format recognitions passed
    //fprintf(stderr, "Error, not recognizing read name formatting: [%s]\n\nIf your data come from SRA, be sure to dump the fastq file like so:\n\n\tSRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra \n\n", name);
    //exit(2);
    
    // now, we'll assume that the user knows what they're doing and just attach the /1 or /2 to the read name.
    
    name[name_len] = '/';
    if (read_type == 1) {
        name[name_len+1] = '1';
    } else {
        name[name_len+1] = '2';
    }
    name[name_len + 2] = '\0';
    
    return;


    
}



static inline void stk_printseq_renamed(const kseq_t *s, int line_len, const char *prefix, int64_t n)
{
	putchar(s->qual.l? '@' : '>');
	if (n >= 0) {
		if (prefix) fputs(prefix, stdout);
		printf("%lld", (long long)n);
	}
    else {
        // trinity mods
        char* name = s->name.s;
        char name_copy [1000];
        strcpy(name_copy, name);
        update_read_name_for_Trinity(name_copy, s->comment.l, s->comment.s);
        fputs(name_copy, stdout);
    }
    
    /*    // trinity mod
	if (s->comment.l) {
		putchar(' '); fputs(s->comment.s, stdout);
	}
    */
	stk_printstr(&s->seq, line_len);
	if (s->qual.l) {
		putchar('+');
		stk_printstr(&s->qual, line_len);
	}
}

static inline void stk_printseq(const kseq_t *s, int line_len)
{
	stk_printseq_renamed(s, line_len, 0, -1);
}

/* 
   64-bit Mersenne Twister pseudorandom number generator. Adapted from:

     http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c

   which was written by Takuji Nishimura and Makoto Matsumoto and released
   under the 3-clause BSD license.
*/

typedef uint64_t krint64_t;

struct _krand_t;
typedef struct _krand_t krand_t;

#define KR_NN 312
#define KR_MM 156
#define KR_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define KR_LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct _krand_t {
	int mti;
	krint64_t mt[KR_NN];
};

static void kr_srand0(krint64_t seed, krand_t *kr)
{
	kr->mt[0] = seed;
	for (kr->mti = 1; kr->mti < KR_NN; ++kr->mti) 
		kr->mt[kr->mti] = 6364136223846793005ULL * (kr->mt[kr->mti - 1] ^ (kr->mt[kr->mti - 1] >> 62)) + kr->mti;
}

krand_t *kr_srand(krint64_t seed)
{
	krand_t *kr;
	kr = malloc(sizeof(krand_t));
	kr_srand0(seed, kr);
	return kr;
}

krint64_t kr_rand(krand_t *kr)
{
	krint64_t x;
	static const krint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
    if (kr->mti >= KR_NN) {
		int i;
		if (kr->mti == KR_NN + 1) kr_srand0(5489ULL, kr);
        for (i = 0; i < KR_NN - KR_MM; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + KR_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        for (; i < KR_NN - 1; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + (KR_MM - KR_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        x = (kr->mt[KR_NN - 1] & KR_UM) | (kr->mt[0] & KR_LM);
        kr->mt[KR_NN - 1] = kr->mt[KR_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
        kr->mti = 0;
    }
    x = kr->mt[kr->mti++];
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
}

#define kr_drand(_kr) ((kr_rand(_kr) >> 11) * (1.0/9007199254740992.0))


int stk_seq(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int c, qual_thres = 0, flag = 0, qual_shift = 33, mask_chr = 0, min_len = 0, max_q = 255;
	unsigned i, line_len = 0;
	int64_t n_seqs = 0;
	double frac = 1.;
	khash_t(reg) *h = 0;
	krand_t *kr = 0;
    

	while ((c = getopt(argc, argv, "N12q:l:Q:aACrn:s:f:M:L:cVUX:SR:")) >= 0) {
		switch (c) {
			case 'a':
			case 'A': flag |= 1; break;
			case 'C': flag |= 2; break;
			case 'r': flag |= 4; break;
			case 'c': flag |= 8; break;
			case '1': flag |= 16; break;
			case '2': flag |= 32; break;
			case 'V': flag |= 64; break;
			case 'N': flag |= 128; break;
			case 'U': flag |= 256; break;
			case 'S': flag |= 512; break;
			case 'M': h = stk_reg_read(optarg); break;
			case 'n': mask_chr = *optarg; break;
			case 'Q': qual_shift = atoi(optarg); break;
			case 'q': qual_thres = atoi(optarg); break;
			case 'X': max_q = atoi(optarg); break;
			case 'l': line_len = atoi(optarg); break;
			case 'L': min_len = atoi(optarg); break;
			case 's': kr = kr_srand(atol(optarg)); break;
			case 'f': frac = atof(optarg); break;
            case 'R': read_type = atoi(optarg); break;
		}
	}
	if (kr == 0) kr = kr_srand(11);
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk seq [options] <in.fq>|<in.fa>\n\n");
		fprintf(stderr, "Options: -q INT    mask bases with quality lower than INT [0]\n");
		fprintf(stderr, "         -X INT    mask bases with quality higher than INT [255]\n");
		fprintf(stderr, "         -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]\n");
		fprintf(stderr, "         -l INT    number of residues per line; 0 for 2^32-1 [%d]\n", line_len);
		fprintf(stderr, "         -Q INT    quality shift: ASCII-INT gives base quality [%d]\n", qual_shift);
		fprintf(stderr, "         -s INT    random seed (effective with -f) [11]\n");
		fprintf(stderr, "         -f FLOAT  sample FLOAT fraction of sequences [1]\n");
		fprintf(stderr, "         -M FILE   mask regions in BED or name list FILE [null]\n");
		fprintf(stderr, "         -L INT    drop sequences with length shorter than INT [0]\n");
		fprintf(stderr, "         -c        mask complement region (effective with -M)\n");
		fprintf(stderr, "         -r        reverse complement\n");
		fprintf(stderr, "         -A        force FASTA output (discard quality)\n");
		fprintf(stderr, "         -C        drop comments at the header lines\n");
		fprintf(stderr, "         -N        drop sequences containing ambiguous bases\n");
		fprintf(stderr, "         -1        output the 2n-1 reads only\n");
		fprintf(stderr, "         -2        output the 2n reads only\n");
		fprintf(stderr, "         -V        shift quality by '(-Q) - 33'\n");
		fprintf(stderr, "         -U        convert all bases to uppercases\n");
		fprintf(stderr, "         -S        strip of white spaces in sequences\n");
        fprintf(stderr, "         -R        read_type 1 (left) or 2 (right).  ie. -R 1 or -R 2\n");
        fprintf(stderr, "\n");
		free(kr);
		return 1;
	}
    
    if (read_type < 1 || read_type > 2) {
        fprintf(stderr, "Error, must specify read type via -R as 1 or 2   ");
        exit(2);
    }

    char* filename = optind < argc && strcmp(argv[optind], "-")  ? argv[optind] : "-";
    
	if (line_len == 0) line_len = UINT_MAX;
	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
	qual_thres += qual_shift;
    
    int reading_read_ret = 1;
	while (reading_read_ret >= 0) {
        

        reading_read_ret = kseq_read(seq);
        
        if (reading_read_ret < -1) {
            // error encountered.  see kseq_read() method
            fprintf(stderr, "Error encountered just after sequence entry[%li]: %s, quals and seq lines dont match in length:\n\n... corrupt file?", n_seqs+1, seq->name.s);
            exit(3);
        }
        else if (reading_read_ret == -1) {
            // done reading file.
            break;
        }
        else if (reading_read_ret == 0) {
            fprintf(stderr, "Error encountered at sequence entry[%li] ... corrupt file?", n_seqs);
            exit(4);
        }

        ++n_seqs;
        //fprintf(stderr, "-processing seq: %li\n", n_seqs);
		
        
        if (seq->seq.l < min_len) continue; // NB: length filter before taking random
		if (frac < 1. && kr_drand(kr) >= frac) continue;
		if (flag & 48) { // then choose odd/even reads only
			if ((flag&16) && (n_seqs&1) == 0) continue;
			if ((flag&32) && (n_seqs&1) == 1) continue;
		}
		if (flag & 512) { // option -S: squeeze out white spaces
			int k;
			if (seq->qual.l) {
				for (i = k = 0; i < seq->seq.l; ++i)
					if (!isspace(seq->seq.s[i]))
						seq->qual.s[k++] = seq->qual.s[i];
				seq->qual.l = k;
			}
			for (i = k = 0; i < seq->seq.l; ++i)
				if (!isspace(seq->seq.s[i]))
					seq->seq.s[k++] = seq->seq.s[i];
			seq->seq.l = k;
		}
		if (seq->qual.l && qual_thres > qual_shift) {
			if (mask_chr) {
				for (i = 0; i < seq->seq.l; ++i)
					if (seq->qual.s[i] < qual_thres || seq->qual.s[i] > max_q)
						seq->seq.s[i] = mask_chr;
			} else {
				for (i = 0; i < seq->seq.l; ++i)
					if (seq->qual.s[i] < qual_thres || seq->qual.s[i] > max_q)
						seq->seq.s[i] = tolower(seq->seq.s[i]);
			}
		}
		if (flag & 256) // option -U: convert to uppercases
			for (i = 0; i < seq->seq.l; ++i)
				seq->seq.s[i] = toupper(seq->seq.s[i]);
		if (flag & 1) seq->qual.l = 0; // option -a: fastq -> fasta
		if (flag & 2) seq->comment.l = 0; // option -C: drop fasta/q comments

		if (flag & 4) { // option -r: reverse complement
			int c0, c1;
			for (i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
				c0 = comp_tab[(int)seq->seq.s[i]];
				c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
				seq->seq.s[i] = c1;
				seq->seq.s[seq->seq.l - 1 - i] = c0;
			}
			if (seq->seq.l & 1) // complement the remaining base
				seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];
			if (seq->qual.l) {
				for (i = 0; i < seq->seq.l>>1; ++i) // reverse quality
					c0 = seq->qual.s[i], seq->qual.s[i] = seq->qual.s[seq->qual.l - 1 - i], seq->qual.s[seq->qual.l - 1 - i] = c0;
			}
		}
		if ((flag & 64) && seq->qual.l && qual_shift != 33)
			for (i = 0; i < seq->qual.l; ++i)
				seq->qual.s[i] -= qual_shift - 33;
		if (flag & 128) { // option -N: drop sequences containing ambiguous bases - Note: this is the last step!
			for (i = 0; i < seq->seq.l; ++i)
				if (seq_nt16to4_table[seq_nt16_table[(int)seq->seq.s[i]]] > 3) break;
			if (i < seq->seq.l) continue;
		}
		stk_printseq(seq, line_len);
	}

    if (n_seqs <= 0) {
        fprintf(stderr, "Error, no records were correctly parsed from %s", filename);
        exit(5);
    }
    
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	free(kr);
	return 0;
}


/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   seqtk <command> <arguments>\n");
	fprintf(stderr, "Version: 1.2-r95-dirty\n\n");
	fprintf(stderr, "Command: seq       common transformation of FASTA/Q\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();

    if (strcmp(argv[1], "seq") == 0) stk_seq(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
