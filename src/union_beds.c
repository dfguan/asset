/*
 * =====================================================================================
 *
 *       Filename:  union_beds.c
 *
 *    Description:  merge two beds file
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


void cord_ctg_push(cord_ctg_t *cc, int n)
{
	if (n >= cc->n) {
		if (n >= cc->m) {
			cc->m = n? n << 1 : 2;
			cc->cord_ctgs = realloc(cc->cord_ctgs, cc->m * sizeof(cord_t));
		}
		cc->cord_ctgs[n] = (cord_t){0,0,0};
		cc->n++;	
	}
}

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
	if (a->s < b->s) 
		return -1;
	else if (a->s == b->s) {
		if (a->e < b->e) 
			return -1;
		else if (a->e == b->e) 
			return 0;
		else 
			return 1;
	} else 
		return 1;
}

cord_ctg_t *col_b1(char *bed_fn, sdict_t *ctgs, bed_hdr_t *hdr)
{
	bed_file_t *bfp = bed_open(bed_fn);
	if (!bfp) return NULL;
	
	cord_ctg_t *cc = calloc(1, sizeof(cord_ctg_t));
	bed_hdr_read(bfp, hdr);
	fprintf(stderr, "test:%s", hdr->desc);
	bed_rec_t r;
	while (bed_read(bfp, &r) >= 0) {
		/*fprintf(stderr, "%s\n", r.ctgn);*/
		int ind = sd_put(ctgs, r.ctgn, 10);
		cord_ctg_push(cc, ind);
		/*fprintf(stderr, "%d\n", ind);*/
		cors tmp = (cors){r.s, r.e}; 
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
		int ind = sd_put(ctgs, r.ctgn, 10); //not check ing here 
		cord_ctg_push(cc, ind);
		cors tmp = (cors){r.s, r.e}; 
		cord_push(&cc->cord_ctgs[ind], &tmp);	
	}	
	bed_close(bfp);
	return 0;
}

cord_ctg_t *union_bs(cord_ctg_t *c1, uint32_t n, sdict_t *ctgs)
{
	cord_ctg_t *cn = calloc(1, sizeof(cord_ctg_t));
	cn->cord_ctgs = calloc(n, sizeof(cord_t));
	
	uint32_t i;
	for (i = 0; i < n; ++i) {
		size_t n_c = c1->cord_ctgs[i].n;
		/*fprintf(stderr, "%u\n", n_c);*/
		cors *vec = c1->cord_ctgs[i].coords;
		uint32_t j;
		/*for ( j = 1; j < n_c; ++j) {*/
			/*fprintf(stderr, "%u\t%u\n", vec[j].s, vec[j].e);*/
		/*}*/
		qsort(vec, n_c, sizeof(cors), cmp);	
		cors a = (cors) {vec[0].s, vec[0].e};
		/*fprintf(stderr, "%u\n", n_c);*/
		/*for ( j = 0; j < n_c; ++j) {*/
			/*fprintf(stderr, "%s\t%u\t%u\n", ctgs->seq[i].name, vec[j].s, vec[j].e);*/
		/*}*/
		for ( j = 1; j <= n_c; ++j) {
			if (j == n_c || vec[j].s > a.e) {
				/*fprintf(stderr, "%s\t%u\t%u\n", ctgs->seq[i].name, vec[j].s, vec[j].e);*/
				cord_push(&cn->cord_ctgs[i], &a);
				if (j < n_c) a = (cors) {vec[j].s, vec[j].e};
			} else {
				a.e = max(vec[j].e, a.e);
			} 
		}
	}
	return cn;
}

int gen_bed(cord_ctg_t *cn, sdict_t *ctgs, bed_hdr_t *hdr)
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
		cors *coords = cn->cord_ctgs[i].coords;
		uint32_t n_coords = cn->cord_ctgs[i].n;
	/*fprintf(stderr, "%u\n", n_coords);*/
		uint32_t j;
		for ( j = 0; j < n_coords; ++j) {
			fprintf(stdout, "%s\t%u\t%u\n", ctgn,coords[j].s, coords[j].e);
		}	
	}
	return 0;
}


int main(int argc, char *argv[])
{
	
	if (argc < 2) {
		fprintf(stderr,"[E::%s] require two bed files each time!\n", __func__); 
		return -1;
	}
	char *b1_fn = argv[1];
	char *b2_fn = argv[2];

	sdict_t *ctgs = sd_init();
	bed_hdr_t h = (bed_hdr_t){0,0};
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collect information from bed file 1\n", __func__);
#endif
	
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
	cord_ctg_t *cn =union_bs(cc, ctgs->n_seq, ctgs);
	cord_ctg_destroy(cc);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] write bed information\n", __func__);
#endif
	gen_bed(cn, ctgs,&h);
	cord_ctg_destroy(cn);	
	sd_destroy(ctgs);
	return 0;
}





