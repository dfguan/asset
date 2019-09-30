/*
 * =====================================================================================
 * *       Filename:  col_hic_lnks.c
 *
 *    Description:  collect links from bam files 
 *
 *        Version:  1.0
 *        Created:  21/10/2018 10:04:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>


#include "sdict.h"
#include "cdict.h"
#include "bamlite.h"
#include "bed.h"
#include "kvec.h"

typedef struct {
	/*int mq:15, rev:1, as:16;*/
	uint32_t s, ns;
   	uint32_t tid:31, rev:1;
}aln_inf_t;

typedef struct {
	uint64_t c1ns:63, qrev:1;
	uint64_t c2ns:63, trev:1;
}hit_t;
typedef struct {
	uint32_t n, m;
	hit_t *ary;
}hit_ary_t;

void hit_ary_push(hit_ary_t *l, hit_t *z)
{
	uint32_t max = -1;
	if (l->n >= l->m) {
		if (l->m > (max >> 1)) {
			fprintf(stderr, "Too many values here\n");
			exit(1);
		} else 
			l->m = l->m ? l->m << 1 : 16;
		l->ary = realloc(l->ary, sizeof(hit_t) * l->m);//beware of overflow
	}
	l->ary[l->n++] = *z;
}

int cmp_hits(const void *a, const void *b)
{
	hit_t *m = (hit_t *)a;
	hit_t *n = (hit_t *)b; //too many branches
	if (m->c1ns > n->c1ns) return 1;	
	else if (m->c1ns < n->c1ns) return -1;
	else if (m->c2ns > n->c2ns) return 1;	
	else if (m->c2ns == n->c2ns) return 0;
	else return -1;

}

uint32_t check_left_half(uint32_t le, uint32_t rs, uint32_t p) // 1 for left half 0 for right half 2 for middle
{
	if (p > le && p < rs) return 2;
	else if (p <= le) return 1;
	else return 0;	
}
/*
graph_t *gen_graph(cdict_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	uint32_t n = ctgs->n_seq << 1;
	uint32_t i, j;
	for ( i = 0; i < n; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			uint32_t ind = sd_get(ctgs, name2) << 1 | c->cnts[j].is_l;
			uint32_t k;
			uint8_t hand_shaking = 0;
			uint32_t ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			if (hand_shaking) fprintf(stderr, "I hand shaking\n");
			if (hand_shaking) add_edge2(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	
		}		
	}	
	

	return g;
}
*/
/*
int get_edge_from_txt(char *links_fn, cdict_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {

		uint32_t ind1 = sd_get(ctgs, r.ctgn);		
		line_n += 1;
		cd_add2(&cds[ind1<<1|r.is_l], r.ctgn2, r.is_l2, r.wt, r.llen);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}
*/


void out_matrix(cdict_t *cds, sdict_t *ctgs, uint32_t n, char *out_fn)
{
	uint32_t i;
	cdict_t *c;

	FILE *fout = out_fn ? fopen(out_fn, "w") : stdout;
	for ( i = 0; i < n; ++i) {
		c = cds + i;
		uint32_t j;
		for ( j = 0; j < c->n_cnt; ++j) {
			if (c->cnts[j].cnt) fprintf(fout, "%s\t%c\t%s\t%c\t%u\t%u\t%u\n", ctgs->seq[i>>1].name, i&1?'+':'-', c->cnts[j].name, j&1?'+':'-', c->cnts[j].cnt, c->cnts[j].snp_n, i & 1 ? ctgs->seq[i>>1].l_snp_n : ctgs->seq[i>>1].r_snp_n);				
		}	
	}
	if (out_fn) fclose(fout);
}

/*  
int core(char *snps_fn, char *edge_fn)
{
	sdict_t *ctgs = init_snps(snps_fn);	
	if (!ctgs) return 1;
	
	fprintf(stderr, "%u\n", ctgs->n_seq);
	cdict_t* cds = calloc(ctgs->n_seq<<1, sizeof(cdict_t)); 
	uint32_t n_cds = ctgs->n_seq<<1;
	
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
	get_edge_txt(edge_fn, cds, ctgs);
	anothernorm(cds, ctgs);
	return 0;
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	cd_set_lim(cds, n_cds); 
	graph_t *g = gen_graph(cds, ctgs);
	process_graph(g);
	
	fprintf(stderr, "enter\n");

	for (i = 0; i < n_cds; ++i) {
		fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);
		cd_destroy(cds +i);	

	} 
	fprintf(stderr, "leave\n");
	if (cds) free(cds);
	graph_destroy(g);
	return 0;

}
*/
int col_enzcuts(hit_ary_t *hit_ary, sdict_t *sd)
{
	size_t i, j;
	
    uint32_t *idxs = calloc(sd->n_seq + 1, sizeof(uint32_t)); //bug when genome larger than 32G
    for (i = 0; i < sd->n_seq; ++i) 
        idxs[i] = (sd->seq[i].len + 7) / 8;
    uint32_t sum = 0, tmp;
    for (i = 0; i <= sd->n_seq; ++i) 
       tmp = idxs[i], idxs[i] = sum, sum += tmp; 
    
    for (i = 0; i <= sd->n_seq; ++i) 
        fprintf(stderr, "%u\n", idxs[i]); 
    hit_t *hs = hit_ary->ary;
	size_t n = hit_ary->n;
    
    uint8_t sets[8] = {1, 2, 4, 8, 16, 32, 64, 128};  
    uint8_t *enzyme_hits = calloc(sum, sizeof(uint8_t));

	for (i = 0, j = 1; j <= n; ++j) {
		if (j == n || hs[i].c1ns != hs[j].c1ns || hs[i].c2ns != hs[j].c2ns || hs[i].qrev != hs[j].qrev || hs[i].trev != hs[j].trev) {
			uint32_t ind1 = hs[i].c1ns >> 32; 
			uint32_t ind2 = hs[i].c2ns >> 32; 
			uint32_t a0s = (uint32_t) hs[i].c1ns; 
			uint32_t a1s = (uint32_t) hs[i].c2ns; 
	        enzyme_hits[idxs[ind1] + (a0s >> 3)] |= sets[a0s & 0x7];	
	        enzyme_hits[idxs[ind2] + (a1s >> 3)] |= sets[a1s & 0x7];	
			i = j;	
		}
	}
    uint8_t *count_bits_tab = malloc(sizeof(uint8_t) * 256);
    int num_to_bits[16] =  {0, 1, 1, 2, 1, 2, 2, 3,
                            1, 2, 2, 3, 2, 3, 3, 4};
    for (i = 0; i < 256; ++i) 
        count_bits_tab[i] = num_to_bits[i & 0xF] + num_to_bits[(i>>4) & 0xF]; 
    
    for ( i = 0; i < sd->n_seq; ++i) {
        sd->seq[i].l_snp_n = sd->seq[i].r_snp_n = 0;       
        uint32_t halfl = sd->seq[i].len >> 1;      
        for ( j = 0; j < ((halfl + 7) >> 3); ++j) 
            sd->seq[i].l_snp_n += count_bits_tab[enzyme_hits[idxs[i] + j]];
            sd->seq[i].l_snp_n += count_bits_tab[enzyme_hits[idxs[i] + j] >> (8 - (halfl & 0x7))]; 
        sd->seq[i].r_snp_n += count_bits_tab[enzyme_hits[idxs[i] + j] & (halfl & 0x7)]; 
        for ( ++j; j < idxs[i+1] - idxs[i]; ++j) 
            sd->seq[i].r_snp_n += count_bits_tab[enzyme_hits[idxs[i] + j]];
    }
    free(count_bits_tab);
    free(enzyme_hits);
    free(idxs);
}


int col_contacts(hit_ary_t *hit_ary, sdict_t *sd, cdict_t *cs)
{
	size_t i, j;
	sdict_t *use_sd = sd;
	hit_t *hs = hit_ary->ary;
	size_t n = hit_ary->n;
	for (i = 0, j = 1; j <= n; ++j) {
		if (j == n || hs[i].c1ns != hs[j].c1ns || hs[i].c2ns != hs[j].c2ns || hs[i].qrev != hs[j].qrev || hs[i].trev != hs[j].trev) {
			
			uint32_t ind1 = hs[i].c1ns >> 32; 
			uint32_t ind2 = hs[i].c2ns >> 32; 
			uint32_t a0s = (uint32_t) hs[i].c1ns; 
			uint32_t a1s = (uint32_t) hs[i].c2ns; 
			uint32_t is_l1 = check_left_half(use_sd->seq[ind1].le, use_sd->seq[ind1].rs, a0s);
			if (is_l1 > 1) return 1; //middle won't be added
			uint32_t is_l2 = check_left_half(use_sd->seq[ind2].le, use_sd->seq[ind2].rs, a1s);
			if (is_l2 > 1) return 1; //middle won't be added
			fprintf(stderr, "%s\t%u\t%s\t%u\n", sd->seq[ind1].name, a0s, sd->seq[ind2].name, a1s);
			cd_add(&cs[ind1<<1|is_l1], use_sd->seq[ind2].name, is_l2, is_l2?use_sd->seq[ind2].l_snp_n:use_sd->seq[ind2].r_snp_n);		
			cd_add(&cs[ind2<<1|is_l2], use_sd->seq[ind1].name, is_l1, is_l1?use_sd->seq[ind1].l_snp_n:use_sd->seq[ind1].r_snp_n);		
			
			i = j;	
		}
	}
}



int col_hits(aln_inf_t *a, int a_cnt, aln_inf_t *f, int f_cnt, sdict_t *ctgs, sdict_t *scfs, hit_ary_t *hit_ary)
{
	if (scfs->n_seq) {
		if (a_cnt == 2) {
			/*fprintf(stderr, "%u\t%u\n", a[0].tid, a[1].tid);*/
			sd_seq_t *sq1 = &ctgs->seq[a[0].tid];
			sd_seq_t *sq2 = &ctgs->seq[a[1].tid];

			uint32_t ind1 = sq1->le; //maybe not well paired up
			uint32_t ind2 = sq2->le;
			if (ind1 == ind2) return 1;
			/*fprintf(stderr, "%s\t%s\t%u\t%u\n", sq1->name, sq2->name, ind1, ind2);*/
			/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
			uint32_t a0s = sq1->l_snp_n == a[0].rev ? sq1->rs + a[0].s : sq1->rs + sq1->len - a[0].s; 
			uint32_t a1s = sq2->l_snp_n == a[1].rev ? sq2->rs + a[1].s : sq2->rs + sq2->len - a[1].s; 
			
			if (ind1 < ind2) {
				uint64_t c1ns = (uint64_t)ind1 << 32 | a0s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind2 << 32 | a1s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, a[0].rev, c2ns, a[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			} else {
				uint64_t c1ns = (uint64_t)ind2 << 32 | a1s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind1 << 32 | a0s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, a[0].rev, c2ns, a[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			}
			return 0;
		} else if (f_cnt == 2){
			sd_seq_t *sq1 = &ctgs->seq[f[0].tid];
			sd_seq_t *sq2 = &ctgs->seq[f[1].tid];
			uint32_t ind1 = sq1->le; //maybe not well paired up
			uint32_t ind2 = sq2->le;
			if (ind1 == ind2) return 1;
			/*fprintf(stderr, "%u\t%u\n", ind1, ind2);*/
			/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
			uint32_t f0s = sq1->l_snp_n == f[0].rev ? sq1->rs + f[0].s : sq1->rs + sq1->len - f[0].s; 
			uint32_t f1s = sq1->l_snp_n == f[1].rev ? sq2->rs + f[1].s : sq2->rs + sq2->len - f[1].s; 
			if (ind1 < ind2) {
				uint64_t c1ns = (uint64_t)ind1 << 32 | f0s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind2 << 32 | f1s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, f[0].rev, c2ns, f[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			} else {
				uint64_t c1ns = (uint64_t)ind2 << 32 | f1s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind1 << 32 | f0s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, f[0].rev, c2ns, f[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			}
			return 0;	
		}
	} else {
		if (a_cnt == 2) {
			uint32_t ind1 = a[0].tid; //maybe not well paired up
			uint32_t ind2 = a[1].tid;
			if (ind1 == ind2) return 1;
			/*fprintf(stderr, "%u\t%u\n", ind1, ind2);*/
			/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
			uint32_t is_l1 = check_left_half(ctgs->seq[ind1].le, ctgs->seq[ind1].rs, a[0].s);
			if (is_l1 > 1) return 1; //middle won't be added
			uint32_t is_l2 = check_left_half(ctgs->seq[ind2].le, ctgs->seq[ind2].rs, a[1].s);
			if (is_l2 > 1) return 1; //middle won't be added
			
			if (ind1 < ind2) {
				uint64_t c1ns = (uint64_t)ind1 << 32 | a[0].s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind2 << 32 | a[1].s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, a[0].rev, c2ns, a[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			} else {
				uint64_t c1ns = (uint64_t)ind2 << 32 | a[1].s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind1 << 32 | a[0].s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, a[0].rev, c2ns, a[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			}
			return 0;
		} else if (f_cnt == 2){
			uint32_t ind1 = f[0].tid;
			uint32_t ind2 = f[1].tid;
			if (ind1 == ind2) return 1;
			/*fprintf(stderr, "%s\t%s\n", r->ctgn1, r->ctgn2)	;*/
			if (ind1 < ind2) {
				uint64_t c1ns = (uint64_t)ind1 << 32 | f[0].s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind2 << 32 | f[1].s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, f[0].rev, c2ns, f[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			} else {
				uint64_t c1ns = (uint64_t)ind2 << 32 | f[1].s; //don't think there will be 2G contig, if happends might be a bug 
				uint64_t c2ns = (uint64_t)ind1 << 32 | f[0].s; //don't think there will be 2G contig, if happends might be a bug 
				hit_t h = (hit_t) {c1ns, f[0].rev, c2ns, f[1].rev}; 
				hit_ary_push(hit_ary, &h);	
			}
			return 0;	
		}
	}
	return 1;
}

int proc_bam(char *bam_fn, int min_mq, sdict_t *ctgs, sdict_t *scfs, hit_ary_t *ha)
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
	/*int i;*/
	/*for ( i = 0; i < h->n_targets; ++i) {*/
		/*char *name = h->target_name[i];*/
		/*uint32_t len = h->target_len[i];*/
		/*if (ws > len) ws = len;*/
		/*uint32_t le = (len - ws) >> 1;*/
		/*uint32_t rs = (len + ws) >> 1;*/
		/*uint32_t lenl, lenr;*/
		/*lenl = lenr = (len - ws) >> 1;*/
		/*sd_put2(ctgs, name, len, le, rs, lenl, lenr);*/
	/*}*/
	/*uint32_t cur_ws;*/
	/*for ( i = 0; i < h->n_targets; ++i) {*/
		/*char *name = h->target_name[i];*/
		/*uint32_t len = h->target_len[i];*/
		/*cur_ws = ws;*/
		/*if (len < (cur_ws << 1)) cur_ws = len >> 1;*/
		/*uint32_t le = cur_ws;*/
		/*uint32_t rs = len - cur_ws + 1;*/
		/*uint32_t lenl, lenr;*/
		/*lenl = lenr = cur_ws;*/
		/*sd_put2(ctgs, name, len, le, rs, lenl, lenr);*/
	/*}*/
	/*if (!ns->ct) { //not initiate yet*/
		/*init_gaps(gap_fn, ns, ctgs, max_ins_len);*/
	/*}*/

	char *cur_qn = 0;
	long bam_cnt = 0;
	int is_set = 0;
	/*aln_inf_t aln[2];*/
	/*int aln_cnt;*/
	
	kvec_t(aln_inf_t) all;
	kv_init(all);
	kvec_t(aln_inf_t) five;
	kv_init(five);

	uint8_t rev;
	uint64_t rdp_counter  = 0;
	uint64_t used_rdp_counter = 0;
	/*fprintf(stderr, "Proc Bam %d\n", __LINE__);*/
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				if (!col_hits(all.a, all.n, five.a, five.n, ctgs, scfs, ha)) ++used_rdp_counter;
				/*aln_cnt = 0;	*/
				/*rev = 0;*/
				/*is_set = 0;*/
				kv_reset(all);
				kv_reset(five);
				if (cur_qn) ++rdp_counter, free(cur_qn); 
				cur_qn = strdup(bam1_qname(b));
			}
			if (b->core.flag & 0x4 || b->core.qual < min_mq) continue; //not aligned
			aln_inf_t tmp;
			tmp.rev = !!(b->core.flag & 0x10);
			/*tmp.nrev = !!(b->core.flag & 0x20);*/
			//only collects five prime
			tmp.tid = b->core.tid;
			/*tmp.ntid = b->core.mtid;*/
			tmp.s = b->core.pos + 1;
			/*tmp.ns = b->core.mpos + 1;*/
			kv_push(aln_inf_t, all, tmp);
			
			uint32_t *cigar = bam1_cigar(b);
			if ((rev && bam_cigar_op(cigar[0]) == BAM_CMATCH) || (!rev && bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CMATCH)) 
				kv_push(aln_inf_t, five, tmp);
			
			/*aln_cnt = (aln_cnt + 1 ) & 1;*/
			/*if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); */
		} else {
			if (!col_hits(all.a, all.n, five.a, five.n, ctgs, scfs, ha)) ++used_rdp_counter;
			if (cur_qn) ++rdp_counter, free(cur_qn); 
			break;	
		}
	}
	fprintf(stderr, "[M::%s] finish processing %lld read pairs %lld (%.2f) passed\n", __func__, rdp_counter, used_rdp_counter, (double)used_rdp_counter/rdp_counter); 
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	kv_destroy(all);
	kv_destroy(five);
	return 0;
}


int chl_col_ctgs(char *bam_fn, sdict_t *ctgs, uint32_t ws)
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
	
	/*ctg_pos_t *d = ctg_pos_init();*/
	/*for ( i = 0; i < h->n_targets; ++i) */
		/*sd_put(ctgs, h->target_name[i], h->target_len[i]);*/
		/*ctg_pos_push(d, i);*/
	//50k
	/*uint32_t cur_ws;*/
	/*for ( i = 0; i < h->n_targets; ++i) {*/
		/*char *name = h->target_name[i];*/
		/*uint32_t len = h->target_len[i];*/
		/*cur_ws = ws;*/
		/*if (len < (cur_ws << 1)) cur_ws = len >> 1;*/
		/*uint32_t le = cur_ws;*/
		/*uint32_t rs = len - cur_ws + 1;*/
		/*uint32_t lenl, lenr;*/
		/*lenl = lenr = cur_ws;*/
		/*sd_put2(ctgs, name, len, le, rs, lenl, lenr);*/
	/*}*/
	int i;
	for ( i = 0; i < h->n_targets; ++i) {
		char *name = h->target_name[i];
		uint32_t len = h->target_len[i];
		uint32_t le = len >> 1;
		uint32_t rs = (len >> 1) + 1;
		uint32_t lenl, lenr;
		lenl = lenr = len >> 1;
		sd_put2(ctgs, name, len, le, rs, lenl, lenr);
	}
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	return 0;
}



/*int aa_10x_hic(char *bam_fn, int min_as, int min_mq, int min_cov, float min_cov_rat, int max_cov, float max_cov_rat)*/
/*int aa_hic(char *bam_fn, int min_as, int min_mq, int min_cov, int max_cov, uint32_t max_ins_len)*/
int col_hic_lnks(char *sat_fn, char **bam_fn, int n_bam, int min_mq, uint32_t win_s, char *out_fn)
{

	/*uint32_t n_cds = ctgs->n_seq<<1;*/
	/*cdict_t* cds = calloc(ctgs->n_seq<<1, sizeof(cdict_t)); */

	sdict_t *ctgs = sd_init();	
	sdict_t *scfs = sd_init();

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] initiate contigs\n", __func__);
#endif
	chl_col_ctgs(bam_fn[0], ctgs, win_s);	
	if (!ctgs) {
		fprintf(stderr, "[E::%s] fail to collect contigs\n", __func__);	
		return 1;
	} 

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif
	
	hit_ary_t *hit_ary = calloc(1, sizeof(hit_ary_t));
	int i;	
	for ( i = 0; i < n_bam; ++i) {
		if (proc_bam(bam_fn[i], min_mq, ctgs, scfs, hit_ary)) {
			return 1;	
		}	
	}
	//sort hit_ary
	if (!hit_ary->ary) {
		fprintf(stderr, "[W::%s] no qualified hits found in the alignments\n", __func__);
		return 1;
	} 
	qsort(hit_ary->ary, hit_ary->n, sizeof(hit_t), cmp_hits);	
	//col joints
    col_enzcuts(hit_ary, ctgs);
	fprintf(stderr, "collect thing");
	sdict_t *_sd = scfs->n_seq ? scfs : ctgs;
	cdict_t *cds = calloc(_sd->n_seq << 1, sizeof(cdict_t));
	for ( i = 0; i < _sd->n_seq << 1; ++i) cd_init(&cds[i]); //be careful with the access way
	col_contacts(hit_ary, _sd, cds);
	
	fprintf(stderr, "finish thing");
	free(hit_ary->ary); free(hit_ary);

	uint32_t n_cds = _sd->n_seq << 1;
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
	out_matrix(cds, _sd, n_cds, out_fn);
	for (i = 0; i < n_cds; ++i)  cd_destroy(cds +i);	
	if (cds) free(cds);
	sd_destroy(ctgs);
	sd_destroy(scfs);
	return 0;

}

int main(int argc, char *argv[])
{
	int c;
	int  min_mq = 10;
	/*uint32_t max_ins_len = 10000;*/
	/*int max_cov = 100, min_cov = 0, min_mq = 0;*/
	/*int min_as = 0;*/
	/*uint32_t max_ins_len = 10000;*/
	uint32_t win_s = 50000;
	char *program;
	char *sat_fn = 0, *out_fn = 0;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "q:w:s:h"))) {
		switch (c) {
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'w':
				win_s = strtol(optarg, NULL, 10);
				break;
			case 's':
				sat_fn  = optarg;
				break;
			case 'o':
				out_fn  = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [options] <BAM_FILE>\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [10]\n");
				fprintf(stderr, "         -w    INT      window size [50000]\n");
				fprintf(stderr, "         -s    STR      sat file\n");
				fprintf(stderr, "         -o    STR      output file\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	char **bam_fn = &argv[optind];
	int n_bam = argc - optind;
	fprintf(stderr, "Program starts\n");	
	col_hic_lnks(sat_fn, bam_fn, n_bam, min_mq, win_s, out_fn);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}

