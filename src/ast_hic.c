/*
 * =====================================================================================
 *
 *       Filename:  aa_hic.c
 *
 *    Description:  assembly assessment with hic 
 *
 *        Version:  1.0
 *        Created:  15/09/2018 15:42:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>


#include "ast.h"
#include "bamlite.h"
#include "sdict.h"
#include "kvec.h"
#include "ksort.h"
#include "bed.h"

typedef struct {
	int mq:15, rev:1, as:15, nrev:1;
	uint32_t s, e, tid;
}aln_inf_t;


typedef struct {
	size_t n, m;
	uint64_t *se;
} se_t;

typedef struct {
	int n;
	se_t *st_ary;
}se_ary_t;

void se_push(se_t *st, uint64_t _se)
{
	uint32_t max = -1;
	if (st->n >= st->m) {
		st->m = st->m ? st->m << 1 : 16;
		st->se = realloc(st->se, sizeof(uint64_t) * st->m);//beware of overflow
	}
	st->se[st->n++] = _se;
}

void sary_destroy(se_ary_t *sary)
{
	if (sary && sary->n) {
		int i;
		for ( i = 0; i < sary->n; ++i ) if (sary->st_ary[i].n) free(sary->st_ary[i].se);	
		free(sary);
	}
	return;
}

#define u64_key(a) ((a))
KRADIX_SORT_INIT(u64, uint64_t, u64_key, 8)

int bin_srch(uint32_t s, uint32_t e, cors *cs, uint32_t nc, uint32_t min_is) //search first greater
{
	cors *c = cs;
	if (!nc) return -1; 
	if (s < c[0].s || s > c[nc-1].e) return -1;
	uint32_t mid;
	uint32_t h = nc, l = 0;
	while (l != h) {
		mid = (l+h)>>1;
		if (c[mid].s <= s) l = mid + 1;
		else 
			h = mid;
	}	
	/*if (e > c[l].e) return 0;*/
	/*e - s - gap < 10kb;*/
	if (e > c[l].e && (l == nc-1 || e < c[l+1].s) && (e - s) < c[l].e - c[l].s + min_is) return 0;
	/*if (e > c[l].e && (l == nc-1 || e < c[l+1].s)) return l;*/
	else
		return -1;
}

int span_gaps(uint32_t s, uint32_t e, cord_t *ct, uint32_t min_is)
{
	cors *cs = ct->coords;
	int nc = ct->n;
	int ind = bin_srch(s, e, cs, nc, min_is);
	if (~ind) return 1;
	else 
	/*if (~ind) {*/
		 /*calculate middle point*/
		/*uint32_t mids, mide; */
		/*if (!ind) {*/
			/*mids = (0 + cs[ind].s) >> 1	;*/
		/*} else */
			/*mids = (cs[ind-1].e + cs[ind].s) >> 1;	*/
		/*if (ind == nc - 1)*/
			/*mide = (cs[ind].e + len) >> 1; */
		/*else*/
			/*mide = (cs[ind].e + cs[ind+1].s) >> 1;*/
		/*if (s >= mids && e <= mide) */
			/*return 1;*/
	/*}  */
	return 0;
}


/*void col_pos(aln_inf_t  *fal, int min_as, int min_mq, uint32_t max_ins_len, ctg_pos_t *d)*/
/*
void col_pos(aln_inf_t  *fal, int min_mq, uint32_t max_ins_len, ns_t *ns, ctg_pos_t *d)
{
	if (fal->mq > min_mq ) {
		uint32_t s = fal->s;
		uint32_t e = fal->e;
		if (e - s < max_ins_len || span_gaps(s, e, &ns->ct[fal->tid])) {
			s = s << 1;
			e = e << 1 | 1; //
			pos_push(&d->ctg_pos[fal->tid], s);
			pos_push(&d->ctg_pos[fal->tid], e); // we init 
		}
	}
}
*/

void col_pos2(aln_inf_t  *fal, ctg_pos_t *d)
{
	/*if (fal->mq > min_mq ) {*/
		uint32_t s = fal->s;
		uint32_t e = fal->e;
		/*fprintf(stderr, "%u\t%u\t%u\t%d\t%u\n", s, e, e - s, min_mq, max_ins_len);*/
		/*if (e - s < max_ins_len || span_gaps(s, e, &ns->ct[fal->tid])) {*/
	/*if (e-s >= max_ins_len) fprintf(stderr, "%u\t%u\t%u\n", s, e, e - s);*/
		s = s << 1;
		e = e << 1 | 1; //
		pos_push(&d->ctg_pos[fal->tid], s);
		pos_push(&d->ctg_pos[fal->tid], e); // we init 
		/*}*/
	/*}*/
}
uint32_t get_target_end(uint32_t *cigar, int n_cigar, uint32_t s)
{
	int i = 0;
	for ( i = 0; i < n_cigar; ++i) {
		uint8_t c  = bam_cigar_opchr(cigar[i]);
		if (c == 'M' || c == 'D') 
			s += cigar[i] >> BAM_CIGAR_SHIFT;
	}	
	return s;
}

void init_gaps(char *gap_fn, ns_t *ns, sdict_t *ctgs)
{
	bed_file_t* bf = bed_open(gap_fn);
	bed_rec_t r;
	ns->ct = calloc(ctgs->n_seq, sizeof(cord_t));
	while (bed_read(bf, &r) >= 0) {
		uint32_t ind = sd_put(ctgs, r.ctgn, 0);
		/*if (r.e - r.s >= min_len) {*/
			cors tmp = (cors){r.s + 1, r.e};
			cord_push(&ns->ct[ind], &tmp);				
		/*}*/
	}
	bed_close(bf);
}

/*  
uint8_t check_alns(aln_inf_t *a, int aln_cnt, uint8_t rev, uint32_t max_is, ns_t *ns)
{
	// 0 pass others alignments with problems
	uint8_t v = 0;
	if (aln_cnt != 2) {
		v = 1;
	} else {
		if (a[0].mq <= 30) v |= (1 << 1);
		if (rev != 1 && rev != 2) v |= (1 << 2);
		if (a[0].e - a[0].s > max_is) {
			v |= (1 << 3);
			if (!span_gaps(a[0].s, a[0].e, &ns->ct[a[0].tid])) v |= (1 << 4); 
		}
	}	
	
	return v;
}
*/

int col_pos(se_ary_t *sary, ctg_pos_t *d) 
{
	int k;
	for ( k = 0; k < sary->n; ++k) {
		size_t i, j; 
		uint64_t *st = sary->st_ary[k].se;
		size_t n = sary->st_ary[k].n;
		for ( i = 0, j = 1; j <= n; ++j) {
			if (j == n || st[i] != st[j]) {
				uint32_t s = st[i] >> 33 << 1;
				uint32_t e = (uint32_t) st[i] | 1;	
				fprintf(stderr, "%u\t%u\n", s, e);
				pos_push(&d->ctg_pos[k], s);		
				pos_push(&d->ctg_pos[k], e);		
				i = j;
			}	
		}		
	}	
	return 0;
}
// this need to be changed to support repeats maybe we should filter out PCR dupliates?

/*int col_suprt(aln_inf_t *all, int all_cnt, aln_inf_t *five, int five_cnt, ns_t *ns, uint32_t max_is, ctg_pos_t *d, sdict_t *ctgs)*/
int col_suprt(aln_inf_t *all, int all_cnt, aln_inf_t *five, int five_cnt, ns_t *ns, uint32_t max_is, se_ary_t *sary, sdict_t *ctgs)
{
	if (all_cnt == 1 ) {
		/*if (all->rev != all->nrev && (all->e - all->s < max_is || span_gaps(all->s, all->e, &ns->ct[all->tid], ctgs->seq[all->tid].len))) {*/
		if (all->e - all->s < max_is || span_gaps(all->s, all->e, &ns->ct[all->tid], max_is)) {
			uint64_t so = (all->s << 1) | (all->rev ? 1 : 0);// allow contig size less than 2G 
			all->e = (all->e << 1) | (all->nrev ? 1 : 0);
			/*fprintf(stderr, "%lx\n", so<<32 | all->e);*/
			se_push(&sary->st_ary[all->tid], (so << 32) | all->e);
			/*uint32_t e = all->e << 1 | 1;*/
			/*pos_push(&d->ctg_pos[all->tid], s);*/
			/*pos_push(&d->ctg_pos[all->tid], e);*/
			return 0;
		}	
	} else if (five_cnt == 1) {
		if (five->e - five->s < max_is || span_gaps(five->s, five->e, &ns->ct[five->tid], max_is)) {
			uint64_t so = (five->s << 1) | (five->rev ? 1 : 0);// allow contig size less than 2G 
			five->e = (five->e << 1) | (five->nrev ? 1 : 0);
			se_push(&sary->st_ary[five->tid], (so << 32) | five->e);
			/*uint32_t s = five->s << 1;*/
			/*uint32_t e = five->e << 1 | 1;*/
			/*pos_push(&d->ctg_pos[five->tid], s);*/
			/*pos_push(&d->ctg_pos[five->tid], e);*/
			return 0;
		}	
	}	
	return 1;
}


/*int proc_bam(char *bam_fn, int min_mq, uint32_t max_ins_len, sdict_t *ctgs, ctg_pos_t *d, char *gap_fn, ns_t *ns)*/
int proc_bam(char *bam_fn, int min_mq, uint32_t max_ins_len, sdict_t *ctgs, se_ary_t *sary, char *gap_fn, ns_t *ns)
{
	
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	fp = bam_open(bam_fn, "r"); //should check if bam is sorted
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open %s\n", __func__, bam_fn);
		return -1;
	}
	
	h = bam_header_read(fp);
	b = bam_init1();
	
	int i;
	for ( i = 0; i < h->n_targets; ++i) {
		sd_put(ctgs, h->target_name[i], h->target_len[i]);
		/*ctg_pos_push(d, i);*/
	}
	if (!sary->n) {
		sary->st_ary = calloc(ctgs->n_seq, sizeof(se_t));
		sary->n = ctgs->n_seq;
	}

	if (!ns->ct) { //not initiate yet
		init_gaps(gap_fn, ns, ctgs);
	}
	


	char *cur_qn = NULL;
       /**cur_bc = NULL;*/
	/*int32_t cur_l = 0;*/
	/*kvec_t(aln_inf_t) fsa = {0, 0, 0}; // fisrt segment*/
	/*kvec_t(aln_inf_t) lsa = {0, 0, 0}; //last segemnt*/
	/*sdict_t* bc_n = sd_init();*/
	/*long bam_cnt = 0;*/
	/*int is_set = 0;*/
	/*aln_inf_t aln[4];*/
	/*int aln_cnt;*/
	/*uint8_t rev, rtv;*/
	
	kvec_t(aln_inf_t) all;
	kvec_t(aln_inf_t) five;
	kv_init(all);
	kv_init(five);
	/*uint64_t err_counter[6] = {0};*/
	uint64_t rdp_counter = 0;
	uint64_t used_rdp_counter = 0;
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				if (!col_suprt(all.a, all.n, five.a, five.n, ns, max_ins_len, sary, ctgs)) ++used_rdp_counter;
				kv_reset(all);
				kv_reset(five);
				if (cur_qn) free(cur_qn); 
				cur_qn = strdup(bam1_qname(b));
			}
			if (b->core.flag & 0x4 || b->core.flag & 0x8 || b->core.qual < min_mq || b->core.isize <= 0) continue; //not aligned
			aln_inf_t tmp; 
			tmp.rev = !!(b->core.flag & 0x10);
			tmp.nrev = !!(b->core.flag & 0x20);
			tmp.s = b->core.pos + 1; //one-based 
			tmp.mq = b->core.qual;
			tmp.tid = b->core.tid;
			tmp.e = tmp.s + b->core.isize - 1;//fully closed	
			kv_push(aln_inf_t, all, tmp);	
			uint32_t *cigar = bam1_cigar(b);
			if ((tmp.rev && bam_cigar_op(cigar[0]) == BAM_CMATCH)|| (!tmp.rev && bam_cigar_op(cigar[b->core.n_cigar - 1]) == BAM_CMATCH))	
			kv_push(aln_inf_t, five, tmp);
				/*} */
			if ((++rdp_counter % 1000000) == 0) fprintf(stderr, "[M::%s] processing %lld read pairs\n", __func__, rdp_counter); 
		} else {
			if (!col_suprt(all.a, all.n, five.a, five.n, ns, max_ins_len, sary, ctgs)) ++used_rdp_counter;
			break;	
		}
	}
	fprintf(stderr, "[M::%s] Total: %lld read pairs, %lld (%.2lf) Passed\n", __func__, rdp_counter, used_rdp_counter, (double)used_rdp_counter/rdp_counter);
/*fprintf(stderr, "Total: %lld read pairs, %lld (%.2lf) Passed %lld (%.2lf) Number Error %lld (%.2lf) MAQ Error, %lld (%.2lf) Direction Error, %lld (%.2lf) IS Error, %lld (%.2lf) GAP Error", rdp_counter, err_counter[0], (double)err_counter[0]/rdp_counter, err_counter[1], (double)err_counter[1]/rdp_counter, err_counter[2], (double)err_counter[2]/rdp_counter, err_counter[3],(double)err_counter[3]/rdp_counter, err_counter[4], (double)err_counter[4]/rdp_counter, err_counter[5], (double) err_counter[5]/rdp_counter);*/
	/*fprintf(stderr, "[M::%s] finish processing %ld bams\n", __func__, bam_cnt); */
	/*while (1) {*/
		//segment were mapped 
		/*if (bam_read1(fp, b) >= 0 ) {*/
			/*fprintf(stderr, "%s\t%s\n",cur_qn, bam1_qname(b));*/
			/*if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {*/
				/*col_pos(fsa.a, fsa.n, lsa.a, lsa.n, min_as, min_mq,max_ins_len, d);*/
				/*if (cur_qn) free(cur_qn); */
				/*cur_qn = strdup(bam1_qname(b));*/
				/*cur_bc = cur_qn + b->core.l_qname - BC_LEN;*/
				/*cur_l = b->core.l_qseq;*/
				/*lsa.n = fsa.n = 0;*/
			/*}*/
			/*if (b->core.flag & 0x4) continue; //not aligned*/
			
			/*aln_inf_t tmp;*/
			/*tmp.rev = !!(b->core.flag & 0x10);*/
			/*uint8_t *s = bam_aux_get(b, "AS");*/
			/*if (s) tmp.as = *(int32_t *)(s+1); else tmp.as = -1;	*/
			/*tmp.s = b->core.pos;*/
			/*tmp.mq = *(uint16_t *)bam1_qual(b);*/
			/*tmp.tid = b->core.tid;*/
			/*uint32_t e = get_target_end(bam1_cigar(b), b->core.n_cigar, tmp.s);			*/
			/*tmp.e = e;	*/
			/*if (tmp.rev) {*/
				/*uint32_t tmp_pos = tmp.s;*/
				/*tmp.s = ctgs->seq[b->core.tid].len - tmp.e;*/
				/*tmp.e = ctgs->seq[b->core.tid].len - tmp_pos;	*/
			/*}*/
			/*if ((b->core.flag & 0x40) && !(b->core.flag & 0x80)) {*/
				/*add reverse or forward infor*/
				/*tmp.tn = (uint32_t)gfa_name2id(g, h->target_name[b->core.tid]);*/
					
				/*uint8_t *s = bam_aux_get(b, "NM");*/
				/*if (s) tmp.nm = *(int32_t *)(s+1); else tmp.nm = -1;*/
				/*kv_push(aln_inf_t,fsa,tmp);				*/
			/*} else if ((b->core.flag & 0x80) && !(b->core.flag & 0x40)) {*/
				/*aln_inf_t tmp;*/
				/*tmp.tn = (uint32_t)gfa_name2id(g, h->target_name[b->core.tid]);*/
				/*tmp.rc = !!(b->core.flag & 0x10);*/
				/*uint8_t *s = bam_aux_get(b, "NM");*/
				/*if (s) tmp.nm = *(int32_t *)(s+1); else tmp.nm = -1;*/
				/*uint8_t *s = bam_aux_get(b, "AS");*/
				/*if (s) tmp.as = *(int32_t *)(s+1); else tmp.as = -1;	*/
				/*tmp.s = b->core.pos;*/
				/*kv_push(aln_inf_t, lsa, tmp); 				*/
			/*}*/
			/*if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); */
		/*} else {*/
			/*col_pos(fsa.a, fsa.n, lsa.a, lsa.n, min_as, min_mq, max_ins_len, d);*/
			/*calc_alns(g, fsa.a, fsa.n, lsa.a, lsa.n, cur_l, fas, sd_put(bc_n,cur_bc, BC_LEN), seg_bc);*/
			/*break;	*/
		/*}*/
	/*}*/
	
	/*kv_destroy(fsa);*/
	/*kv_destroy(lsa);	*/
	kv_destroy(all);
	kv_destroy(five);
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	/*sd_destroy(bc_n);	*/
	return 0;
}

int sort_se(se_ary_t *sary)
{
	int i;
	for ( i = 0; i < sary->n; ++i ) {
		if (sary->st_ary[i].n) {
			uint64_t *st = sary->st_ary[i].se;
			uint64_t *ed = st + sary->st_ary[i].n;
			radix_sort_u64(st, ed);  
		}
	}
	return 0;
}

/*int aa_10x_hic(char *bam_fn, int min_as, int min_mq, int min_cov, float min_cov_rat, int max_cov, float max_cov_rat)*/
/*int aa_hic(char *bam_fn, int min_as, int min_mq, int min_cov, int max_cov, uint32_t max_ins_len)*/
int aa_hic(char **bam_fn, int n_bam, char *gap_fn, int min_mq, int min_cov, int max_cov, uint32_t max_ins_len, char *out_dir)
{
	sdict_t *ctgs = sd_init();
		
	int i;	
	ns_t *ns = calloc(1, sizeof(ns_t));
	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif
	se_ary_t *sary = calloc(1, sizeof(se_ary_t));
	for ( i = 0; i < n_bam; ++i) {
		if (proc_bam(bam_fn[i], min_mq, max_ins_len, ctgs, sary, gap_fn, ns)) {
			return -1;	
		}	
	}	
	
	ns_destroy(ns);	
	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting positions\n", __func__);
#endif
	ctg_pos_t *d = ctg_pos_init();
	sort_se(sary);
	for ( i = 0; i < ctgs->n_seq; ++i ) ctg_pos_push(d, i);
	col_pos(sary, d);
	sary_destroy(sary);

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] calculating coverage\n", __func__);
#endif
	float tmp;
	cov_ary_t *ca = cal_cov(d, ctgs, &tmp);

	ctg_pos_destroy(d);
	if (!ca) {
		fprintf(stderr, "[W::%s] low quality alignments\n", __func__);
		return 0;	
	}
		
	/*sel_sup_reg(ca, min_cov_rat, min_cov, max_cov_rat, max_cov, ctgs);*/
	char *type = "HC";
	char *desc = "hic data";
#ifdef PRINT_COVERAGE
	print_coverage_wig(ca, ctgs, type, 1024, out_dir);
	print_coverage_stat(ca, ctgs, type, out_dir);
	print_base_coverage(ca, ctgs, type, out_dir);
#endif	
	sel_sup_reg(ca, min_cov, max_cov, ctgs, type, desc);
		
	cov_ary_destroy(ca, ctgs->n_seq); //a little bit messy
	sd_destroy(ctgs);

	fprintf(stderr, "Program finished successfully\n");
	return 0;

}

int main(int argc, char *argv[])
{
	int c;
	int max_cov = 10000000, min_cov = 7, min_mq = 0;
	int min_as = 0;
	uint32_t max_ins_len = 15000;
	/*int max_cov = 100, min_cov = 0, min_mq = 0;*/
	/*int min_as = 0;*/
	/*uint32_t max_ins_len = 10000;*/

	char *r;
	char *out_dir = ".";
	while (~(c=getopt(argc, argv, "c:C:q:s:L:O:h"))) {
		switch (c) {
			case 'c':
				min_cov = atoi(optarg); 
				break;
			case 'C':
				max_cov = atoi(optarg); 
				break;
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'L':
				max_ins_len = strtol(optarg, &r, 10);
				break;
			case 'O':
				out_dir = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: aa_hic [options] <GAP_BED> <BAM_FILEs>\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c    INT      minimum coverage [7]\n");
				fprintf(stderr, "         -C    INT      maximum coverage [inf]\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [0]\n");
				/*fprintf(stderr, "         -s    INT      minimum alignment score [0]\n");*/
				fprintf(stderr, "         -L    INT      maximum insertion length, gap excluded [15000]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require at least one bed and bam file!\n", __func__); goto help;
	}
	char *gap_fn = argv[optind++];
	char **bam_fn = &argv[optind];
	int n_bam = argc - optind;
	fprintf(stderr, "Program starts\n");	
	aa_hic(bam_fn, n_bam,  gap_fn, min_mq, min_cov, max_cov, max_ins_len, out_dir);
	return 0;	
}




