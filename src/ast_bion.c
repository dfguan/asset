/*
 * =====================================================================================
 *
 *       Filename:  aa_bionano.c
 *
 *    Description:  assess assembly with bionano data
 *
 *        Version:  1.0
 *        Created:  16/09/2018 18:35:49
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
#include <getopt.h>

#include "ast.h"
#include "amap.h"
#include "kvec.h"
#include "ksort.h"

typedef struct {
	float cov;
	int diff;
	int s, e;
}diff_t;

#define diff_key(a) ((a).diff)
KRADIX_SORT_INIT(dif, diff_t, diff_key, 8)

int cal_bound(diff_t *df, int n_df, int *ub, int *lb)
{
	if (n_df < 2) return -1;//not enough number;
	else if (n_df < 3) {
		*lb = -1;
		*ub = 0x7FFFFFFF;	
	} else {
		radix_sort_dif(df, df + n_df); 
		int q1, q3; //Tukey's hinges
		uint32_t s_2, n; //number for each part
		n = n_df >> 1;
		s_2 = n;
		if (n_df & 1) 
			++n, ++s_2;
		 
		if (n & 1) 
			q1 = df[n >> 1].diff, q3 = df[s_2 + (n>>1)].diff;
		else 
			q1 = (df[(n >> 1) - 1].diff + df[n>>1].diff ) >> 1, q3 = (df[(n >> 1) - 1 + s_2].diff + df[(n >> 1) + s_2].diff) >> 1;
		*ub = (5*q3 - 3*q1) >> 1;
		*lb = 5*q1 > 3*q3 ? (5*q1 - 3*q3) >> 1 : 0;
	}
	return 0;
}

sdict_t *key_to_ctg_name(char *key_fn)
{	
	sdict_t *ctgs = sd_init();
	
	amap_file_t *fp = amap_open(key_fn);
		
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open file: %s\n", __func__, key_fn);
		return NULL;	
	}
	key_entry r;
	int n = 0;
	while (amap_read_k(fp, &r) >= 0) {
		//assume key is sorted and no jumping numbers // not sure if it is okay
		if (!n) {
			++n;
			continue;
		}// skip first line
		sd_put(ctgs, r.ctg_name, (uint32_t)r.ctg_len);		
		++n;
	}
	return ctgs;
}



ctg_pos_t *col_r(char *rcmap, int n_ctg)
{
	ctg_pos_t *rctg_p = ctg_pos_init();	
	int i;
	for ( i = 0; i < n_ctg; ++i) ctg_pos_push(rctg_p, i);
	amap_file_t *fp = amap_open(rcmap);
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open file: %s\n", __func__, rcmap);
		return NULL;	
	}
	rqmap_entry r;
		/*fprintf(stderr,"%u\n", n_ctg);	*/
	while (amap_read_rq(fp, &r) >= 0)  {
		/*fprintf(stderr,"%u\n", r.cmap_id);	*/
		pos_push(&rctg_p->ctg_pos[r.cmap_id-1], (uint32_t)r.pos);	// don't worry ctg_pos has been initiated 
	}
	return rctg_p;
}

typedef struct {
	uint32_t pos;
	float cov;
}qmap_inf;

typedef struct {
	qmap_inf *q;
	int n, m;
}qmap_ary_t;

void qmap_ary_push(qmap_ary_t *q, qmap_inf *t)
{
	if (q->n >= q->m) {
		q->m = q->m ? q->m << 1 : 2;
		q->q = realloc(q->q, sizeof(qmap_inf) * q->m);
	}
	q->q[q->n++] = *t; 
}

void release(qmap_ary_t *q, int n)
{
	int i;
	if (q) {
		for ( i = 0; i <= n; ++i) {
			if (q[i].q) free(q[i].q);
		}
		free(q);	
	}
}

/*void *myrealloc(void *old, size_t n, size_t s)*/
/*{*/
	/*void *new = calloc(n, s);*/
	/*memcpy()*/
/*}*/


qmap_ary_t *col_q(char *qcmap, int* n_qry)
{
	int mcnt = 2;
	qmap_ary_t *qa = calloc(mcnt, sizeof(qmap_ary_t));
	
	amap_file_t *fp = amap_open(qcmap);
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open file: %s\n", __func__, qcmap);
		return NULL;	
	}
	rqmap_entry r;
	while (amap_read_rq(fp, &r) >= 0) {
		if (r.cmap_id >= mcnt) {
			int tmp = mcnt;
			mcnt = r.cmap_id << 1;
			qmap_ary_t *new_qa = calloc(mcnt, sizeof(qmap_ary_t));
			memcpy(new_qa, qa, sizeof(qmap_ary_t) * tmp);
			/*release(qa, mcnt);*/
			if (qa) free(qa);
			qa = new_qa;
			/*mctg_id = r.cmap_id;*/
		}	
		qmap_inf tmp_inf = (qmap_inf){(uint32_t)r.pos, r.cov};
		qmap_ary_push(&qa[r.cmap_id], &tmp_inf);	
	} 
	*n_qry = mcnt;
	return qa;
}

typedef struct {
	diff_t *dfs;
	int		n, m;
}diff_ary_t ;//can we package this data structure? 

void diff_ary_push(diff_ary_t *d, diff_t *item)
{
	if (d->n >= d->m) {
		d->m = d->m ? d->m << 1 : 2;
		d->dfs = realloc(d->dfs, sizeof(diff_t) * d->m);
	}
	d->dfs[d->n++] = *item; 
}

typedef struct {
	int s_ind, q_ind;
}cords;
typedef struct {
	cords *c;
	int n, m;
}cords_ary_t; 

void cords_ary_reset(cords_ary_t *c)
{
	c->n = 0;
}

void cords_ary_push(cords_ary_t *c, cords *t)
{
	if (c->n >= c->m) {
		c->m = c->m ? c->m << 1 : 2;
		c->c = realloc(c->c, sizeof(cords) * c->m);
	}
	c->c[c->n++] = *t;	

}

void cords_ary_destroy(cords_ary_t *c)
{
	if (c && c->c) free(c->c), free(c);
}	

void col_cors(char *s, cords_ary_t *c)
{
	cords_ary_reset(c);
	int i; 
	char *q, *p, *r;
	/*int is_open;*/
	int len_s = strlen(s);

	/*fprintf(stderr, "%s\t%d\n", s, len_s);*/
	for ( p = s + 1, i = 1; i < len_s; ++i) {
		if (s[i] == ',') {
				s[i] = 0;
				q = s + i + 1;
		} else if (s[i] == ')') {
			s[i]=0;
			cords tmp = (cords){strtol(p, &r, 10), strtol(q, &r, 10)};			
			cords_ary_push(c, &tmp);
		} else if (s[i] == '(') p = s + i + 1;
	}
}

diff_ary_t *col_diff_cov2(char *xmap, ctg_pos_t *rmap, qmap_ary_t *qmap, sdict_t *ctgs, float min_conf, int n_ctg)
{
	diff_ary_t *diff_cov = calloc(n_ctg, sizeof(diff_ary_t));	
		
	amap_file_t *fp = amap_open(xmap);
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open file: %s\n", __func__, xmap);
		return NULL;	
	}
	xmap_entry r;
	cords_ary_t *c = calloc(1, sizeof(cords_ary_t));
	while (amap_read_x(fp, &r) >= 0)  {
		if (r.conf > min_conf) {
			col_cors(r.alignment, c);
			int i, rs_ind, re_ind, qs_ind, qe_ind;
			/*fprintf(stderr, "%d\n", c->n);*/
			for ( i = 0; i < c->n - 1; ++i) {
				rs_ind =c->c[i].s_ind;
				re_ind = c->c[i+1].s_ind;
				qs_ind = c->c[i].q_ind;
				qe_ind = c->c[i+1].q_ind;
				if (r.ori == '-') {
					int tmp = qs_ind;
					qs_ind = qe_ind;
					qe_ind = tmp;
				}		
				uint32_t s = rmap->ctg_pos[r.ref_id-1].p[re_ind -1], e = rmap->ctg_pos[r.ref_id-1].p[rs_ind - 1];
				uint32_t diff_ref =  s - e;
				/*fprintf(stderr, "qry_id:%d\t%d\t%d\n",r.qry_id, qs_ind, qe_ind);		*/
				uint32_t diff_qry = qmap[r.qry_id].q[qe_ind - 1].pos - qmap[r.qry_id].q[qs_ind - 1].pos;
				int diff = diff_ref > diff_qry ? diff_ref - diff_qry : diff_qry - diff_ref;
				int ndiff = (float)diff / diff_ref * ctgs->seq[r.ref_id - 1].len;
				float cv = 0;
				int j;
				for ( j = qs_ind; j < qe_ind; ++j) cv += qmap[r.qry_id].q[j-1].cov;
				cv = cv / (qe_ind - qs_ind);
				/*fprintf(stderr, "qry_id:%d\t%d\t%d\tref_id:%d\t%d\t%d\t%d\n",r.qry_id, qs_ind, qe_ind, r.ref_id, rs_ind, re_ind, diff);	*/
				/*fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", ctgs->seq[r.ref_id-1].name, s, e, r.qry_id, qs_ind, qe_ind, diff_qry);	*/
				fprintf(stderr, "%s\t%d\t%d\t%d\n", ctgs->seq[r.ref_id-1].name, e, s,diff);	
				diff_t diff_tmp = (diff_t) {cv, ndiff, rs_ind - 1, re_ind - 1};
				diff_ary_push(&diff_cov[r.ref_id - 1], &diff_tmp);
			}
		}
	}
	cords_ary_destroy(c);
	return diff_cov;
}

diff_ary_t *col_diff_cov(char *xmap, ctg_pos_t *rmap, qmap_ary_t *qmap, float min_conf, int n_ctg)
{
	diff_ary_t *diff_cov = calloc(n_ctg, sizeof(diff_ary_t));	
		
	amap_file_t *fp = amap_open(xmap);
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open file: %s\n", __func__, xmap);
		return NULL;	
	}
	xmap_entry r;
	cords_ary_t *c = calloc(1, sizeof(cords_ary_t));
	while (amap_read_x(fp, &r) >= 0)  {
		if (r.conf > min_conf) {
			col_cors(r.alignment, c);
			int i, rs_ind, re_ind, qs_ind, qe_ind;
			/*fprintf(stderr, "%d\n", c->n);*/
			for ( i = 0; i < c->n - 1; ++i) {
				rs_ind =c->c[i].s_ind;
				re_ind = c->c[i+1].s_ind;
				qs_ind = c->c[i].q_ind;
				qe_ind = c->c[i+1].q_ind;
				if (r.ori == '-') {
					int tmp = qs_ind;
					qs_ind = qe_ind;
					qe_ind = tmp;
				}		
			
				uint32_t s = rmap->ctg_pos[r.ref_id-1].p[re_ind -1], e = rmap->ctg_pos[r.ref_id-1].p[rs_ind - 1];
				uint32_t diff_ref =  s - e;
				/*fprintf(stderr, "qry_id:%d\t%d\t%d\n",r.qry_id, qs_ind, qe_ind);		*/
				uint32_t diff_qry = qmap[r.qry_id].q[qe_ind - 1].pos - qmap[r.qry_id].q[qs_ind - 1].pos;
				int diff = diff_ref > diff_qry ? diff_ref - diff_qry : diff_qry - diff_ref;
				float cv = 0;
				int j;
				for ( j = qs_ind; j < qe_ind; ++j) cv += qmap[r.qry_id].q[j-1].cov;
				cv = cv / (qe_ind - qs_ind);
				/*fprintf(stderr, "qry_id:%d\t%d\t%d\tref_id:%d\t%d\t%d\t%d\n",r.qry_id, qs_ind, qe_ind, r.ref_id, rs_ind, re_ind, diff);	*/
				diff_t diff_tmp = (diff_t) {cv, diff, rs_ind - 1, re_ind - 1};
				diff_ary_push(&diff_cov[r.ref_id - 1], &diff_tmp);
			}
		}
	}
	cords_ary_destroy(c);
	return diff_cov;
}



void sel_sup_reg_bionano(diff_ary_t *diff_cov, sdict_t *ctgs, ctg_pos_t *rmap, float min_cov)
{
	int i;
	int lb, ub;
	uint8_t *rblk = NULL;
	int  n_rblk = 0;
	fprintf(stdout, "track name=\"%s\" description=\"%s\"\n", "BN", "Bionano data");	
	for ( i = 0; i < ctgs->n_seq; ++i) {
		int j;  
		for ( j = 0; j < diff_cov[i].n; ++j) { // too many leavel arrays //change later
		   int u, v;
		   u = diff_cov[i].dfs[j].s;
		   v = diff_cov[i].dfs[j].e;
		   /*fprintf(stderr, "%s\t%u\t%u\t%d\n", ctgs->seq[i].name, rmap->ctg_pos[i].p[u], rmap->ctg_pos[i].p[v], diff_cov[i].dfs[j].diff);*/
		}
		if (cal_bound(diff_cov[i].dfs, diff_cov[i].n, &ub, &lb)) continue;
		/*uint32_t k; for (k = 0; k < diff_cov[i].n; ++k) fprintf(stderr, "%s\t%d\n",ctgs->seq[i].name, diff_cov[i].dfs[k].diff);*/
		/*fprintf(stderr,"enter\n");*/
		int n_ctg_pos = rmap->ctg_pos[i].n;
		if (n_ctg_pos > n_rblk) {
			if(rblk) free(rblk);
			n_rblk = rmap->ctg_pos[i].n;
			rblk = calloc(n_rblk, sizeof(uint8_t));
		} else 
			memset(rblk, 0, sizeof(char)* n_ctg_pos);
		for ( j = 0; j < diff_cov[i].n; ++j) { // too many leavel arrays //change later  
			if (diff_cov[i].dfs[j].cov > min_cov && diff_cov[i].dfs[j].diff > lb && diff_cov[i].dfs[j].diff < ub) {
			int k;
			/*fprintf(stderr, "%s\t%u\t%u\t%d\n", ctgs->seq[i].name, diff_cov[i].dfs[j].s, diff_cov[i].dfs[j].e, diff_cov[i].dfs[j].diff);*/
			for (k = diff_cov[i].dfs[j].s; k < diff_cov[i].dfs[j].e; ++k) 
				rblk[k] = 1;	
			}
		}
		for ( j = 0; j < n_ctg_pos - 1; ++j) {
			int z;
			if(rblk[j])	{
				for (z = j; z < n_ctg_pos - 1 && rblk[z] == 1; ++z);	
				fprintf(stdout, "%s\t%u\t%u\n", ctgs->seq[i].name, rmap->ctg_pos[i].p[j] - 1, rmap->ctg_pos[i].p[z]);
				j = z;
			} 
		}
		/*fprintf(stderr,"leave\n");*/
	} 
	/*fprintf(stderr, "finishing\n");*/
}

int aa_bionano(char *xmap_fn, char *rmap_fn, char *qmap_fn, char *key_fn, float min_conf, float min_cov) 
{
	fprintf(stderr, "Program starts\n");	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contig names and contig length\n", __func__);
#endif
	sdict_t *ctgs = key_to_ctg_name(key_fn);
	if (!ctgs)
		return -1;
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting reference cmap information\n", __func__);
#endif
	ctg_pos_t *rctg_p = col_r(rmap_fn, ctgs->n_seq);
	if (!rctg_p) {
		fprintf(stderr, "[E::%s] no reference information\n", __func__);	
		return -1;
	} 
	int n_qry;
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting query cmap information\n", __func__);
#endif
	qmap_ary_t *qctg_p = col_q(qmap_fn, &n_qry);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting distance differences for each contig\n", __func__);
#endif
	diff_ary_t *dif_cov = col_diff_cov2(xmap_fn, rctg_p, qctg_p, ctgs, min_conf, ctgs->n_seq);
	release(qctg_p, qctg_p->n);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] select supported regions\n", __func__);
#endif

	sel_sup_reg_bionano(dif_cov, ctgs, rctg_p, min_cov);	
	fprintf(stderr, "Program finished successfully\n");
	return 0;
}

int main(int argc, char *argv[])
{
	int c;
	float  min_cov = 3, min_conf = 0;
	while (~(c=getopt(argc, argv, "c:s:h"))) {
		switch (c) {
			case 'c':
				min_cov = atof(optarg); 
				break;
			case 's':
				min_conf = atof(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: aa_bion [options] <REF_CMAP> <QUERY_CMAP> <XMAP> <KEY_FN>\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c    FLOAT      minimum coverage [3.0]\n");
				fprintf(stderr, "         -s    FLOAT      minimum alignment confidence [0.0]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 4 > argc) {
		fprintf(stderr,"[E::%s] require query and reference cmap file, alignment xmap file and key file mapping contig names to contig ID!\n", __func__); goto help;
	}
	char *rmap_fn = argv[optind++];
	char *qmap_fn = argv[optind++];
	char *xmap_fn = argv[optind++];
	char *key_fn = argv[optind];

	aa_bionano(xmap_fn, rmap_fn, qmap_fn, key_fn, min_conf, min_cov);
	
	return 0;	
}
