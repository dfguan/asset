//adapted from hengli
#ifndef EG_BSEQ_H
#define EG_BSEQ_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

struct eg_bseq_file_s;
typedef struct eg_bseq_file_s eg_bseq_file_t;

typedef struct {
	int l_seq, rid;
	char *name, *seq, *qual;
} eg_bseq1_t;

eg_bseq_file_t *eg_bseq_open(const char *fn);
void eg_bseq_close(eg_bseq_file_t *fp);
eg_bseq1_t *eg_bseq_read2(eg_bseq_file_t *fp, int chunk_size, int with_qual, int frag_mode, int *n_);
eg_bseq1_t *eg_bseq_read(eg_bseq_file_t *fp, int chunk_size, int with_qual, int *n_);
eg_bseq1_t *eg_bseq_read_frag(int n_fp, eg_bseq_file_t **fp, int chunk_size, int with_qual, int *n_);
int eg_bseq_destroy(eg_bseq1_t *egb, int n);
int eg_bseq_eof(eg_bseq_file_t *fp);

extern unsigned char seq_nt4_table[256];
extern unsigned char seq_comp_table[256];

static inline int eg_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int eg_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = eg_qname_len(s1);
	l2 = eg_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

static inline void eg_revcomp_bseq(eg_bseq1_t *s)
{
	int i, t, l = s->l_seq;
	for (i = 0; i < l>>1; ++i) {
		t = s->seq[l - i - 1];
		s->seq[l - i - 1] = seq_comp_table[(uint8_t)s->seq[i]];
		s->seq[i] = seq_comp_table[t];
	}
	if (l&1) s->seq[l>>1] = seq_comp_table[(uint8_t)s->seq[l>>1]];
	if (s->qual)
		for (i = 0; i < l>>1; ++i)
			t = s->qual[l - i - 1], s->qual[l - i - 1] = s->qual[i], s->qual[i] = t;
}

#ifdef __cplusplus
}
#endif

#endif
