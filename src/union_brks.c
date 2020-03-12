/*
 * =====================================================================================
 *
 *       Filename:  union_brks.c
 *
 *    Description:  merge contig and scaffold breaks
 *
 *        Version:  1.0
 *        Created:  26/09/2018 15:29:15
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
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "ast.h"
#include "bed.h"
/*#include "kseq.h"*/
/*KSTREAM_INIT(gzFile, gzread, 0x10000)*/


typedef struct {
	cord_t *cord_ctgs;
	size_t n,m; 
} cord_ctg_t;



void cord_ctg_destroy(cord_ctg_t *cc)
{
	if (cc) {
		size_t i;
		for (i = 0; i < cc->n; ++i) {
			if(cc->cord_ctgs[i].coords) free(cc->cord_ctgs[i].coords);	
		}
		free(cc);	
	}

}


int cmp(const void *p, const void *q)
{
	cors *a = (cors *)p;
	cors *b = (cors *)q;
	/*int r = (a->s < b->s || (a->s == b->s && a->e < b->e));  //don't return r here only got 1 or 0;*/
	uint32_t as = a->s >> 1;
	uint32_t bs = b->s >> 1;
	if (as < bs) 
		return -1;
	else if (as == bs) {
		if (a->e < b->e) 
			return -1;
		else if (a->e == b->e) 
			return 0;
		else 
			return 1;
	} else 
		return 1;
}

uint32_t gap_sz(uint32_t s, uint32_t e, cord_t *ct) //search first greater
{
	cors *c = ct->coords;
	uint32_t nc = ct->n;
	uint32_t gap_size = 0;
	if (!nc) return gap_size; 
	if (s < c[0].s || s > c[nc-1].e) return gap_size;
	uint32_t mid;
	uint32_t h = nc, l = 0;
	while (l != h) {
		mid = (l+h)>>1;
		if (c[mid].s <= s) l = mid + 1;
		else 
			h = mid;
	}	
	for (;l < nc && c[l].e <= e; ++l) gap_size += (c[l].e - c[l].s + 1);  
	fprintf(stderr, "gapsz: %d\t%d\t%u\n", s, e, gap_size);
	return gap_size;
}
//near but not span the gap 
int near_gap(uint32_t s, uint32_t e, cors *cs, uint32_t nc, uint32_t allow_gaps) 
{
	cors *c = cs;
	if (!nc) return -1; 
	/*if (s < c[0].s || s > c[nc-1].e) return -1;*/
	uint32_t mid;
	uint32_t h = nc, l = 0;
	while (l != h) {
		mid = (l+h)>>1;
		if (c[mid].s < s) l = mid + 1;
		else 
			h = mid;
	}	
	/*if (e > c[l].e) return 0;*/
	/*e - s - gap < 10kb;*/
	fprintf(stderr, "nearg l: %d %d %d %d %d\n", l, s, e, c[l].s, c[l].e);
	/*if (l == 0) return 0;*/
	if (l == 0) {
		if (e < allow_gaps) return 1;
		else return 0;
	}
	if (e <= c[l].s && (s + allow_gaps > c[l].s || c[l-1].e + allow_gaps > e)) return 1;
	/*if (e > c[l].e && (l == nc-1 || e < c[l+1].s)) return l;*/
	else
		return 0;
}
// the last bit of s is used for gap
cord_ctg_t *col_b1(char *bed_fn, sdict_t *ctgs, bed_hdr_t *hdr)
{
	bed_file_t *bfp = bed_open(bed_fn);
	if (!bfp) return NULL;
	
	cord_ctg_t *cc = calloc(1, sizeof(cord_ctg_t));
		
	cc->cord_ctgs = calloc(ctgs->n_seq, sizeof(cord_t));

	bed_hdr_read(bfp, hdr);
	bed_rec_t r;
	while (bed_read(bfp, &r) >= 0) {
		/*fprintf(stderr, "%s\n", r.ctgn);*/
		int ind = sd_get(ctgs, r.ctgn);
		/*fprintf(stderr, "%d\n", ind);*/
		cors tmp = (cors){r.s << 1, r.e}; 
		cord_push(&cc->cord_ctgs[ind], &tmp);	
	}	
	bed_close(bfp);
	return cc;	
}


int col_b2(char *bed_fn, sdict_t *ctgs, cord_ctg_t *cc)
{
	bed_file_t *bfp = bed_open(bed_fn);
	if (!bfp) return -1;
	
	bed_rec_t r;
	while (bed_read(bfp, &r) >= 0) {
		int ind = sd_get(ctgs, r.ctgn); //not check ing here 
		cors tmp = (cors){r.s << 1 | 1, r.e}; 
		cord_push(&cc->cord_ctgs[ind], &tmp);	
	}	
	bed_close(bfp);
	return 0;
}

cord_ctg_t *union_bs(cord_ctg_t *c1, uint32_t n, sdict_t *ctgs, uint32_t alwgap, ns_t *nt)
{
	cord_ctg_t *cn = calloc(1, sizeof(cord_ctg_t));
	cn->cord_ctgs = calloc(n, sizeof(cord_t));
	
	uint32_t i; uint32_t isscf;
	for (i = 0; i < n; ++i) {
		size_t n_c = c1->cord_ctgs[i].n;
		if (!n_c) continue;
		/*fprintf(stderr, "%u\n", n_c);*/
		cors *vec = c1->cord_ctgs[i].coords;
		uint32_t j = 0;
		/*for ( j = 1; j < n_c; ++j) {*/
			/*fprintf(stderr, "%u\t%u\n", vec[j].s, vec[j].e);*/
		/*}*/
		cord_t *gapct = &nt->ct[i];
		qsort(vec, n_c, sizeof(cors), cmp);	
		cors a = (cors) {vec[0].s >> 1, vec[0].e};
		isscf = vec[0].s & 0x1;
		/*fprintf(stderr, "%u\n", n_c);*/
		/*for ( j = 0; j < n_c; ++j) {*/
			/*fprintf(stderr, "%s\t%u\t%u\n", ctgs->seq[i].name, vec[j].s, vec[j].e);*/
		/*}*/
		fprintf(stderr, "%d: %s\t%u\t%u\n", j, ctgs->seq[i].name, vec[j].s>>1, vec[j].e);
		for ( j = 1; j <= n_c; ++j) {
			fprintf(stderr, "%d: %s\t%u\t%u\n", j, ctgs->seq[i].name, vec[j].s>>1, vec[j].e);
			if (j == n_c || (vec[j].s >> 1) > a.e + alwgap || (!(vec[j].s & 0x1) && gap_sz(a.e - 1, vec[j].s >> 1, gapct))) {
				a.s = a.s << 1 | isscf;
				fprintf(stderr,"insert: %d\t%d\n", a.s >> 1, a.e);
				cord_push(&cn->cord_ctgs[i], &a);
				if (j < n_c) a = (cors) {vec[j].s>>1, vec[j].e}, isscf = vec[j].s & 0x1;
			} else {
				a.e = max(vec[j].e, a.e);
				if (!isscf) isscf = vec[j].s & 0x1;
			} 
		}
	}
	return cn;
}

int gen_bed(cord_ctg_t *cn, sdict_t *ctgs, bed_hdr_t *hdr, int brks_sz, ns_t *ns)
{
	//print header
	if (!hdr->type) hdr->type = strdup("UN");
	if (!hdr->desc) hdr->desc = strdup("unkown data type");
	fprintf(stdout, "track name=\"%s\" description=\"%s\"\n", hdr->type, hdr->desc);	
	//print body 
	int n = ctgs->n_seq;
	int i;
	for ( i = 0; i < n; ++i) {
		char *ctgn = ctgs->seq[i].name;
        uint32_t ctglen = ctgs->seq[i].len;
		cors *coords = cn->cord_ctgs[i].coords;
		uint32_t n_coords = cn->cord_ctgs[i].n;
	/*fprintf(stderr, "%u\n", n_coords);*/
		uint32_t j;
		for ( j = 0; j < n_coords; ++j) {
            if (!(brks_sz && !(coords[j].s & 0x1) && (near_gap(coords[j].s >> 1, coords[j].e, ns->ct[i].coords, ns->ct[i].n, brks_sz)))) 
                fprintf(stdout, "%s\t%u\t%u\t%c\n", ctgn,coords[j].s >> 1, coords[j].e, coords[j].s&0x1? 'S':'C');
		}	
	}
	return 0;
}

void init_gaps(char *gap_fn, ns_t *ns, sdict_t *ctgs, uint32_t min_gapsz)
{
	bed_file_t* bf = bed_open(gap_fn);
	bed_rec_t r;
	while (bed_read(bf, &r) >= 0) {
		uint32_t ind = sd_put(ctgs, r.ctgn, r.e);
		ctgs->seq[ind].len = r.e;
		if (r.e - r.s >= min_gapsz) {
			ns_push(ns, ind);
			cors tmp = (cors){r.s, r.e};
			cord_push(&ns->ct[ind], &tmp);				
		}
	}
	bed_close(bf);
}

int main(int argc, char *argv[])
{
	
    int c;    
	int option = 0; //the way to calculate molecule length //internal parameters not allowed to adjust by users
    uint32_t allow_gaps = 10000;
    int brks_sz = 1000;
	int min_gapsz = 0;
	char *scaf_brks = 0, *ctg_brks = 0;
	while (~(c=getopt(argc, argv, "x:g:s:c:h"))) {
		switch (c) {
			case 'g':
				allow_gaps = atoi(optarg); 
				break;
			case 's':
				min_gapsz = atoi(optarg); 
				break;
			case 'x':
				brks_sz = atoi(optarg); 
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: union_bs [options] <GAP.BED> <CONTIG-LEVEL-BREAKS> <SCAFFOLD-LEVEL-BREAKS>\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -g          gap size for chaining two breaks [10000]\n");
				fprintf(stderr, "         -x          exclude breaks near N kb to the ends of contigs [1000]\n");
				fprintf(stderr, "         -s          minimum gap size to sepearate two contigs in a scaffold [0]\n");
				fprintf(stderr, "         -h          help\n");
				return 1;	
		}		
	}
	if (optind + 3 > argc) {
		fprintf(stderr,"[E::%s] require a gap file and two bed files !\n", __func__); 
        goto help;
	}
	char *gapfn = argv[optind++];
	char *b1_fn = argv[optind++];
	char *b2_fn = argv[optind];

	sdict_t *ctgs = sd_init();
	bed_hdr_t h = (bed_hdr_t){0,0};
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collect information from bed file 1\n", __func__);
#endif
	ns_t *ns = calloc(1, sizeof(ns_t));
	init_gaps(gapfn, ns, ctgs, min_gapsz);

	cord_ctg_t *cc = col_b1(b1_fn, ctgs, &h);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collect information from bed file 2\n", __func__);
#endif
	if (col_b2(b2_fn, ctgs, cc)) {
		return -1;		
	}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] combine information from bed file 2\n", __func__);
#endif
	cord_ctg_t *cn =union_bs(cc, ctgs->n_seq, ctgs, allow_gaps, ns);
	cord_ctg_destroy(cc);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] write bed information\n", __func__);
#endif
	gen_bed(cn, ctgs,&h, brks_sz, ns);
	cord_ctg_destroy(cn);	
	sd_destroy(ctgs);
	ns_destroy(ns);
	return 0;
}





