/*
 * =====================================================================================
 *
 *       Filename:  gen_punchlist.c
 *
 *    Description:  generate a punch list
 *
 *        Version:  1.0
 *        Created:  26/09/2018 18:48:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "bed.h"
#include "sdict.h"
#include "ast.h"

typedef struct {	
	uint32_t s:31, isgap:1, e;
	uint32_t n_tech:16, tech:16;
}bed_inf_t;

typedef struct {
	bed_inf_t *bi;
	size_t n,m;
	char *ctgn;
	bed_inf_t last_one;
	char *last_ctgn;
	int has_last_one;
} bed_ary_t;

void bed_ary_destroy(bed_ary_t *ba)
{
	if (ba) {
		if (ba->bi) free(ba->bi);
		if (ba->ctgn) free(ba->ctgn);
		if (ba->last_ctgn) free(ba->last_ctgn);
		free(ba);
	}
}

int check_near_gap(uint32_t s, uint32_t e, cors *cs, uint32_t nc)
{
	cors *c = cs;
	if (!nc) return 0; 
	if (s < cs[1].s) {
        if (e < cs[1].s) return 0; // not in or cross a gap
        else return 1;
    }
	if (s > cs[nc - 2].e) return 0;
	uint32_t mid;
	uint32_t h = nc - 1, l = 1;
	while (l != h) {
		mid = (l+h)>>1;
		if (c[mid].e < s) l = mid + 1;
		else 
			h = mid;
	}	
    if (e < c[mid].s) return 0;
    else return 1;
}

uint32_t bin_srch(uint32_t s, cors *cs, uint32_t nc) //search first greater
{
	cors *c = cs;
	if (!nc) return -1; 
	if (s < cs[1].s) return 1;
	if (s > cs[nc - 2].s) return nc-1;
	uint32_t mid;
	uint32_t h = nc - 1, l = 1;
	while (l != h) {
		mid = (l+h)>>1;
		if (c[mid].s <= s) l = mid + 1;
		else 
			h = mid;
	}	
	return l;
	/*if (e > c[l].e) return 0;*/
	/*e - s - gap < 10kb;*/
	/*if (e > c[l].e && (l == nc-1 || e < c[l+1].s) && (e - s) < c[l].e - c[l].s + 10000) return 0;*/
	/*if (e > c[l].e && (l == nc-1 || e < c[l+1].s)) return l;*/
	/*else*/
		/*return -1;*/
}

int read_batch_bed(bed_file_t *bfp, bed_ary_t *t)
{
	bed_rec_t r;
	t->n = 0;
	if (t->has_last_one) {
		if (t->n >= t->m) {
			t->m = t->m ? t->m << 1 : 2;
			t->bi = realloc(t->bi, sizeof(bed_inf_t) *t->m);
		}
		t->bi[t->n++] = t->last_one;
		if (t->ctgn) free(t->ctgn);
		t->ctgn = strdup(t->last_ctgn);	
	}
	uint16_t tec;
	int ret;
	while((ret = bed_read(bfp, &r)) >= 0) {
		if (!t->ctgn) { //first time running
			t->ctgn = strdup(r.ctgn);
			if (r.tech) 
				tec = *(uint16_t *)r.tech; 
			else 
				tec = 0;
			bed_inf_t tmp = (bed_inf_t) {r.s, r.isgap, r.e, r.n_tec, tec};			   
			if (t->n >= t->m) {
				t->m = t->m ? t->m << 1 : 2;
				t->bi = realloc(t->bi, sizeof(bed_inf_t) * t->m);
			}
			t->bi[t->n++] = tmp;
		} else if (!strcmp(t->ctgn, r.ctgn)) {
			if (r.tech) 
				tec = *(uint16_t *)r.tech; 
			else 
				tec = 0;
			/*fprintf(stderr, "me: %s\t%s\t%s\n", r.ctgn, t->ctgn, r.tech);*/

			bed_inf_t tmp = (bed_inf_t) {r.s, r.isgap, r.e, r.n_tec, tec};			   
			if (t->n >= t->m) {
				t->m = t->m ? t->m << 1 : 2;
				t->bi = realloc(t->bi, sizeof(bed_inf_t) * t->m);
			}
			t->bi[t->n++] = tmp;
		}else {
			if (r.tech) 
				tec = *(uint16_t *)r.tech; 
			else 
				tec = 0;
			t->last_one = (bed_inf_t) {r.s, r.isgap, r.e, r.n_tec, tec};			   
			if (t->last_ctgn) free(t->last_ctgn);
			t->last_ctgn = strdup(r.ctgn);
			t->has_last_one = 1;
			break;
		}	
	}
	return ret;
}

int check_ends(uint32_t s, uint32_t e, cord_t *ct) //check near the both ends of a contig/chromosomes
{
	uint32_t p = bin_srch(s, ct->coords, ct->n);	
	if (~p) { //not 
		uint32_t i;
		for ( i = p - 1; i < p + 1; ++i) {
			uint32_t cur_s  = ct->coords[p-1].s;
			uint32_t cur_e = ct->coords[p-1].e + 10000;
			cur_s = cur_s > 10000 ? cur_s - 10000 : 0;
			if (s >= cur_s && e <= cur_e) return 0;
		}
	}	
	return 1;

	/*if (e <= 50000 || s + 50000 >= len) return 1;*/
	/*else */
		/*return 0;*/
}

void init_gaps(char *gap_fn, ns_t *ns, sdict_t *ctgs)
{
	bed_file_t* bf = bed_open(gap_fn);
	bed_rec_t r;
	while (bed_read(bf, &r) >= 0) {
		uint32_t ind = sd_put(ctgs, r.ctgn, 0);
		ctgs->seq[ind].len = r.e;
		/*if (r.e - r.s >= max_ins_len) {*/
		ns_push(ns, ind);
		cors tmp = (cors){r.s, r.e}; //don't change to one based
		cord_push(&ns->ct[ind], &tmp);				
		/*}*/
	}
	bed_close(bf);
}

int wrt_pchlst2(bed_ary_t *bd, sdict_t *ctgs, ns_t *ns, int isctg)
{
	uint32_t s = -1, e = -1;
	size_t n = bd->n;
	size_t i;	
	bed_inf_t *a = bd->bi;	
	char *ctgn = bd->ctgn;
	uint32_t ctgid = sd_get(ctgs, ctgn);
	uint32_t ctglen = ctgs->seq[ctgid].len;
	/*cord_t *ct = &ns->ct[ctgid];*/
    /*int neargap;*/
	/*uint32_t ctglen = a[n-1].e;*/
	for ( i = 0; i < n; ++i) {
		if (!a[i].n_tech) {
			if (a[i].isgap && !isctg) {
				if (a[i].s && a[i].e != ctglen)
				fprintf(stdout, "%s\t%u\t%u\n", ctgn, a[i].s, a[i].e);	
			}
			else if (!a[i].isgap && isctg)	
				fprintf(stdout, "%s\t%u\t%u\n", ctgn, a[i].s, a[i].e);	
		}
	}
	/*if (isctg) {*/
		/*for ( i = 0; i < n; ++i) {*/
			/*if (!a[i].n_tech && !a[i].isgap) {*/
				/*size_t j;*/
				/*for (j = i + 1; j < n; ++j) */
					/*if (a[j].n_tech || a[j].isgap) */
						/*break;*/
				/*fprintf(stdout, "%s\t%u\t%u\n", ctgn, a[i].s, a[j-1].e);	*/
				/*i = j;*/
			/*}*/
		/*} */
	
	/*} else {*/
		/*for ( i = 0; i < n; ++i) {*/
			/*if (!a[i].n_tech && !a[i].isgap) {*/
				/*size_t j;*/
				/*for (j = i + 1; j < n; ++j) */
					/*if (a[j].n_tech || a[j].isgap) */
						/*break;*/
				/*fprintf(stdout, "%s\t%u\t%u\n", ctgn, a[i].s, a[j-1].e);	*/
				/*i = j;*/
			/*}*/
		/*} */
	/*}*/
	/*for ( i = 0; i < n; ++i) {*/
		/*if (!a[i].n_tech) {*/
            /*size_t j;*/
            /*if (a[i].isgap) neargap = 1;*/
            /*else neargap = 0;*/
            /*for (j = i + 1; j < n; ++j) {*/
                /*if (!a[j].n_tech) {*/
                    /*if (a[j].isgap) neargap = 1; */
                /*} else */
                    /*break;*/
            /*}*/
            /*if (isctg) {*/
                /*if (!neargap)  */
					/*fprintf(stdout, "%s\t%u\t%u\n", ctgn, a[i].s, a[j-1].e);	*/
            /*} else {*/
                /*if (neargap) */
					/*fprintf(stdout, "%s\t%u\t%u\n", ctgn, a[i].s, a[j-1].e);	*/
            /*}*/
            /*i = j;*/
		/*}*/
	/*} */
	return 0;
}

int wrt_pchlst(bed_ary_t *bd, sdict_t *ctgs, ns_t *ns)
{
	uint32_t s = -1, e = -1;
	size_t n = bd->n;
	size_t i;	
	bed_inf_t *a = bd->bi;	
	char *ctgn = bd->ctgn;
	uint32_t ctgid = sd_get(ctgs, ctgn);
	cord_t *ct = &ns->ct[ctgid];
	/*uint32_t ctglen = a[n-1].e;*/
	for ( i = 0; i < n; ++i) {
		if (a[i].n_tech < 2) {
			s = i, e = i; //index
			size_t j;
			for (j = i+1; j < n; ++j) {
				if ((a[j].n_tech < 2 && j != n-1 && a[j+1].n_tech <= 2)||(a[j].n_tech == 2 && j != n-1 && a[j+1].n_tech < 2)) {
					e = j;
				} else {
					fprintf(stdout, "%s\t%u\t%u\t", ctgn, a[s].s, a[e].e);	
					size_t z;
					for ( z = s; z <= e; ++z) {
						if (a[z].n_tech != 1) 
							fprintf(stdout,"E%u:%u", a[z].n_tech, a[z].e - a[z].s);
						else {
							uint32_t o = a[z].tech;
							/*fprintf(stdout, "tech:%hx\n",o);*/
							fprintf(stdout, "%s:%u", (char *)&o, a[z].e - a[z].s);
						} 
						if (z != e) fprintf(stdout,",");
					}	
					if (!check_ends(a[s].s, a[e].e, ct)) 
						fprintf(stdout, "\tCE");
					fprintf(stdout,"\n");
					i = j + 1;
					break;
				}
			}
			if (j == n) {
				fprintf(stdout, "%s\t%u\t%u\t", ctgn, a[s].s, a[e].e);	
				size_t z;
				for ( z = s; z <= e; ++z) {
					if (a[z].n_tech != 1) 
						fprintf(stdout,"E%u:%u", a[z].n_tech, a[z].e - a[z].s);
					else {
						uint32_t o = a[z].tech;
						fprintf(stdout, "%s:%u", (char *)&o, a[z].e - a[z].s);
					} 
					if (z != e) fprintf(stdout,",");
				}	
				if (!check_ends(a[s].s, a[e].e, ct)) 
						fprintf(stdout, "\tCE");
				fprintf(stdout,"\n");
				break;
			}
		}
	} 
	return 0;
}

int gen_pchlst2(char *bed_fn, sdict_t *ctgs, ns_t *ns, int isctg)
{
	
	bed_file_t *bfp = bed_open(bed_fn);
	bed_ary_t *ba = calloc(1, sizeof(bed_ary_t));	
		
	while(read_batch_bed(bfp, ba) >= 0) wrt_pchlst2(ba, ctgs, ns, isctg);  
	bed_ary_destroy(ba);
	return 0;
}
int gen_pchlst(char *bed_fn, sdict_t *ctgs, ns_t *ns)
{
	
	bed_file_t *bfp = bed_open(bed_fn);
	bed_ary_t *ba = calloc(1, sizeof(bed_ary_t));	
		
	while(read_batch_bed(bfp, ba) >= 0) wrt_pchlst(ba, ctgs, ns);  
	bed_ary_destroy(ba);
	return 0;
}

int main(int argc, char *argv[])
{
    int isctg = 0;
    int c;    
	int option = 0; //the way to calculate molecule length //internal parameters not allowed to adjust by users
	while (~(c=getopt(argc, argv, "ch"))) {
		switch (c) {
			case 'c':
				isctg = 1; 
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: pchlst [options] <GAP_BED><BED_FILE>\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c          search for contig-level error [False]\n");
				fprintf(stderr, "         -h          help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require a evidence accumulated and gap bed file\n", __func__); goto help;
	}
	char *gap_fn = argv[optind++];
	
    ns_t *ns = calloc(1, sizeof(ns_t));
	sdict_t *ctgs = sd_init();
	init_gaps(gap_fn, ns, ctgs);

    char *bed_fn = argv[optind];
	gen_pchlst2(bed_fn, ctgs, ns, isctg);
	ns_destroy(ns);
	sd_destroy(ctgs);
	return 0;
}


