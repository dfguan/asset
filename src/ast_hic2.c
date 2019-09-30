/*
 * =====================================================================================
 *
 *       Filename:  build_graph.c
 *
 *    Description:  build graph with links information
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:39:30
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
#include <math.h>

#include "bed.h"
#include "sdict.h"
#include "cdict.h"
#include "kvec.h"

sdict_t *col_ctgs(char *fn)
{
	bed_file_t *bf = bed_open(fn);
	if (!bf) return 0;
	sdict_t *ctgs = sd_init();		
	bed_rec_t r;
    //compare string before colon
    uint32_t scf_idx = 0;
    uint32_t ctg_idx = 0;
    char *cur_scf = 0; 
    uint32_t ish = 0;
    //aux 1 tail 2 head 3 head & tail 0 mid
	while (bed_read(bf, &r) >= 0)  {
        //find colon, transfer to 0 
        char *p = strchr(r.ctgn, ':');
        *p = 0; 
        if (!cur_scf || strcmp(cur_scf, r.ctgn) != 0) { // a new one
            ++scf_idx;
            if (cur_scf) free(cur_scf);
            cur_scf = strdup(r.ctgn); 
            ish = 1;
            if (ctg_idx != 0) ctgs->seq[ctg_idx - 1].aux |= 1; //as tail
        } 
        *p = ':';
		sd_put2(ctgs, r.ctgn, r.s, scf_idx, 0, 0, 0);
        if (ish) ish = 0, ctgs->seq[ctg_idx].aux = 2;
        else ctgs->seq[ctg_idx].aux = 0; 
        ++ctg_idx;
    }
    if (ctg_idx != 0) ctgs->seq[ctg_idx - 1].aux |= 1; // as tail 
	bed_close(bf);
	return ctgs;
}
int print_lnks(cdict2_t *cds, sdict_t *ctgs) 
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
        char *name1 = ctgs->seq[i].name;
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
	    fprintf(stderr, "%d\n", c->n_cnt);
		for (j = 0; j < c->n_cnt; ++j) {
			char *name2 = c->cnts[j].name; 
            fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\t%lf\n", name1,  name2, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3], c->cnts[j].ncnt);
		}
	}
    return 0;
}
int get_links(char *links_fn, cdict2_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1;
		if (r.is_l) 
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, 0);		
		else
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, 0, r.llen);		
		if (r.is_l2)
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, r.rlen, 0);		
		else
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, 0, r.rlen);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd2_add(&cds[ind1], !r.is_l, r.ctgn2, !r.is_l2, r.wt);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}


int norm_links(cdict2_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		/*char *name1 = ctgs->seq[i>>1].name;*/
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
		uint32_t icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            c->cnts[j].ncnt = (float) icnt / ctgs->seq[ctg2_idx].len;
            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			/*fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));*/
		}
	}
	return 0;
}

//from  https://github.com/bcgsc/arcs
//
float norm_cdf(int x, float p, int n) {
    float mean = n * p;
    float sd = sqrt(n * p * (1 - p));
    return 0.5 * (1 + erf((x - mean)/(sd * sqrt(2))));
}

int get_max(uint32_t *a, int n)
{
    int i, idx = -1;
    uint32_t max = 0;
    for ( i = 0; i < n; ++i) 
        if (a[i] > max) idx = i, max = a[i]; 
    return idx;
}

int get_fmax(float *a, int n)
{
    int i, idx = -1;
    float max = 0.0;
    for ( i = 0; i < n; ++i) 
        if (a[i] > max) idx = i, max = a[i]; 
    return idx;
}
inline int ismax(uint32_t *a, int idx, int n)
{
    int i;
    for ( i = 0; i < n; ++i) if (a[i] > a[idx]) break;
    return i == n;
}

uint32_t *find_breaks2(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        uint32_t scf_idx = ctgs->seq[i].le;
        char *name = ctgs->seq[i].name;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
            uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                ctgn += 1;
                fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                if (ctgn >= 10) break;
                //check connect to head to tail 
                
                uint32_t tl, hd, tl2, hd2;
                tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];
                hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];
                /*hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];*/
                /*tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];*/
                //is successor?
                if ((tl > hd && norm_cdf(tl, 0.5, tl + hd)> 0.95) || ismax(c->cnts[j].cnt, 2, 4) || ismax(c->cnts[j].cnt, 3, 4)) {
                    if  (!fgt) ++fgt;
                    else continue;
                    /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                    if (ctg_idx2 != ctg_idx + 1) {
                        fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                        fprintf(stderr, "%s:n1\n", name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } 
                
                } else if ((tl < hd && norm_cdf(hd, 0.5, tl + hd) > 0.95) || ismax(c->cnts[j].cnt, 0, 4) || ismax(c->cnts[j].cnt, 1, 4)) {
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                        fprintf(stderr, "%s:n0\n", name);
                        /*fprintf(stderr, "%s\t%s\n", name, ctgs->seq[ctg_idx+1].name);*/
						/*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                        if (!(ctgs->seq[ctg_idx].aux&0x2)) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);                        } 
                
                } else {
                    fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                    if (ctgn >= 10) break;
                if (ctg_idx2 > ctg_idx) {
                    if (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
						fprintf(stderr, "%s:n1\n", name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } else { // validate direction
                        uint32_t tl, hd, tl2, hd2;
                        tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];
                        hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];
                        hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];
                        tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];
                        if (!((tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95) && (hd2 > tl2 || norm_cdf(tl2, 0.5, tl2 + hd2) <= 0.95)))  {
                        /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
							fprintf(stderr, "%s:m1\n", name);
                            kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                        }
                    }
                } else {
                    if (flt == 0) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
						/*fprintf(stderr, "%s:n0\n", name);*/
                        fprintf(stderr, "%s\t%s\n", name, ctgs->seq[ctg_idx+1].name);
						/*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } else {
                        uint32_t tl, hd, tl2, hd2;
                        tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];
                        hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];
                        hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];
                        tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];
                        if (!((hd > tl || norm_cdf(tl, 0.5, tl + hd) <= 0.95) && (tl2 > hd2 || norm_cdf(hd2, 0.5, tl2 + hd2) <= 0.95))) {
                        /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
							fprintf(stderr, "%s:m0\n", name);
							/*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                        }
                    }
                } 
                }
                }
                    /*
                if (tl > hd || ismax(c->cnts[j].cnt,2, 4))  {
                    if  (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
                        fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                        fprintf(stderr, "%s:n1\n", name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } 
                } else{
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                        fprintf(stderr, "%s:n0\n", name);
                        if (!(ctgs->seq[ctg_idx].aux&0x2)) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    }
                }
            }
                */
            if (fgt && flt) break;
		}
        if (fgt == 0 && (ctgs->seq[ctg_idx].aux & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].aux >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        if (ctgs->seq[ctg_idx].aux & 1) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
	}
    *n_brks = v.n;
	return v.a;
}


uint32_t *find_breaks3(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        uint32_t scf_idx = ctgs->seq[i].le;
        char *name = ctgs->seq[i].name;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
            uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                //is_suc 
                fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                ctgn += 1;
                    if (ctgn >= 10) break;
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (ctg_idx2 > ctg_idx) {
                    /*fprintf(stderr, "suc: %s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                    if (!(tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95)) continue; 
                    if (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
						fprintf(stderr, "%s:n1\n",  name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } 
                    /*else { // validate direction*/
                        /*if (!((tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95) && (hd2 > tl2 || norm_cdf(tl2, 0.5, tl2 + hd2) <= 0.95)))  {*/
                        /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
							/*fprintf(stderr, "%s:m1\n", name);*/
                            /*kv_push(uint32_t, v, ctg_idx << 1 | 1);         */
                        /*}*/
                    /*}*/
                } else {
                    /*fprintf(stderr, "pre: %s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                    if (!(tl < hd || norm_cdf(tl, 0.5, tl + hd) <= 0.95)) continue; 
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                        fprintf(stderr, "%s:n0\n", name);
						/*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } 
                    /*else {*/
                        /*if (!((hd > tl || norm_cdf(tl, 0.5, tl + hd) <= 0.95) && (tl2 > hd2 || norm_cdf(hd2, 0.5, tl2 + hd2) <= 0.95))) {*/
                        /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
							/*fprintf(stderr, "%s:m0\n", name);*/
                            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                            /*kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         */
                        /*}*/
                    /*}*/
                }
                /*
                if (ctg_idx2 > ctg_idx) {
                    if (tl < hd) continue; 
                    if (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
						fprintf(stderr, "%s:n1\n", name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } else { // validate direction
                        if (!((tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95) && (hd2 > tl2 || norm_cdf(tl2, 0.5, tl2 + hd2) <= 0.95)))  {
							fprintf(stderr, "%s:m1\n", name);
                            kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                        }
                    }
                } else {
                    if (tl > hd) continue; 
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        fprintf(stderr, "%s:n0\n", name);
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } else {
                        uint32_t tl, hd, tl2, hd2;
                        tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];
                        hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];
                        hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];
                        tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];
                        if (!((hd > tl || norm_cdf(tl, 0.5, tl + hd) <= 0.95) && (tl2 > hd2 || norm_cdf(hd2, 0.5, tl2 + hd2) <= 0.95))) {
							fprintf(stderr, "%s:m0\n", name);
                            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                        }
                    }
                } 
            */
            }
            if (fgt && flt) break;
		}
        if (fgt == 0 && (ctgs->seq[ctg_idx].aux & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].aux >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        if (ctgs->seq[ctg_idx].aux & 1) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
	}
    *n_brks = v.n;
	return v.a;
}

uint32_t *find_breaks4(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        uint32_t scf_idx = ctgs->seq[i].le;
        char *name = ctgs->seq[i].name;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
            uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                //is_suc 
                ctgn += 1;
                /*if (ctgn >= 10) break;*/
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (ctgn <= 30) fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\t%f\t%s\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3], tl > hd?norm_cdf(tl, 0.5, tl + hd): norm_cdf(hd, 0.5, tl + hd), tl > hd?"TAIL":"HEAD");
                if (ctg_idx2 > ctg_idx && ctg_idx2 == ctg_idx + 1) {
                    if (!((tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95))) { 
                        // this is a precessor 
                        ++flt;
                        fprintf(stderr, "%s:%s\tm1\n",  name, ctgs->seq[ctg_idx2].name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } else 
                        ++fgt;
                } else if (ctg_idx2 < ctg_idx && ctg_idx2 + 1 == ctg_idx) {
                    if (!((hd > tl || norm_cdf(tl, 0.5, tl + hd) <= 0.95))) { //successor
                        ++fgt;
                        fprintf(stderr, "%s:%s\tm0\n", ctgs->seq[ctg_idx2].name, name);
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } else 
                        ++flt; //precedssor
                } else {
                    //check if successor or precessor
                    if (hd > tl || !(get_max(c->cnts[j].cnt, 4) & 0x2)) { //this is a precedssor 
                        if (!flt) {
                            if (!(ctgs->seq[ctg_idx].aux & 0x2)) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                            ++flt; 
                        }
                    } else {
                        if (!fgt) {
                            kv_push(uint32_t, v, (ctg_idx) << 1 | 1);         
                            ++fgt; 
                        }
                    } 
                } 
                /*
                if (ctg_idx2 > ctg_idx) {
                    if (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
                        if (flt && (hd > tl || !(get_max(c->cnts[j].cnt, 4) & 0x2))) { --fgt; continue;}
						fprintf(stderr, "%s:%s\tn1\n",  name, ctgs->seq[ctg_idx2].name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } else { // validate direction
                    }
                } else {
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        if (fgt && (tl > hd || (get_max(c->cnts[j].cnt, 4) & 0x2))) { --flt; continue;}
						fprintf(stderr, "%s:%s\tn0\n",  ctgs->seq[ctg_idx2].name, name);
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } else {
                    }
                }
            */
            }
            if (fgt && flt) break;
		}
        if (fgt == 0 && (ctgs->seq[ctg_idx].aux & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].aux >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        if (ctgs->seq[ctg_idx].aux & 1) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
	}
    *n_brks = v.n;
	return v.a;
}


uint32_t *find_breaks6(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        uint32_t scf_idx = ctgs->seq[i].le;
        char *name = ctgs->seq[i].name;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
            uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                //is_suc 
                if (ctg_idx2 < ctg_idx) {} 
            
            
            }
            if (fgt && flt) break;
		}
        if (fgt == 0 && (ctgs->seq[ctg_idx].aux & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].aux >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        if (ctgs->seq[ctg_idx].aux & 1) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
	}
    *n_brks = v.n;
	return v.a;
}



uint32_t *find_breaks7(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        uint32_t scf_idx = ctgs->seq[i].le;
        char *name = ctgs->seq[i].name;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
        float density[4];
        uint32_t susp_hd = -1, susp_tl = -1;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
            uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                //is_suc 
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                ctgn += 1;
                /*for (z = 0; z < 4; ++z) density[z] = (float)c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[ctg_idx].r_snp_n : ctgs->seq[ctg_idx].l_snp_n) / (z&0x1 ? ctgs->seq[ctg_idx2].r_snp_n : ctgs->seq[ctg_idx2].l_snp_n);*/
                for (z = 0; z < 4; ++z) density[z] = (float)c->cnts[j].cnt[z]/((z >> 1 ? ctgs->seq[ctg_idx].r_snp_n : ctgs->seq[ctg_idx].l_snp_n) + (z&0x1 ? ctgs->seq[ctg_idx2].r_snp_n : ctgs->seq[ctg_idx2].l_snp_n));
                if (ctgn <= 50) fprintf(stderr, "%u\t%s\t%u\t%s\t%u\t%u\t%u\t%u\t%f\t%f\t%f\t%f\n", ctg_idx, name, ctg_idx2, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3], density[0], density[1], density[2], density[3]);

                    /*if (ctgn >= 10) break;*/
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                {
                        if (ctg_idx2 > ctg_idx) {
                            /*if (get_max(c->cnts[j].cnt, 4) > 1) {*/
                            /*if (get_fmax(density, 4) > 1) {*/
                            /*if (tl > hd || get_max(c->cnts[j].cnt, 4) > 1) {*/
                            if (tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95) { //a very loose condition for successor otherwise cause many false positives
                                if (fgt) continue;
                                ++fgt;
                                if (ctg_idx2 != ctg_idx + 1) {
                                    fprintf(stderr, "n1: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                    if (norm_cdf(tl, 0.5, tl + hd) > 0.95) kv_push(uint32_t, v, ctg_idx << 1 | 1);     
                                    else susp_tl = j;
                                } 
                            } // a successor
                            else {
                                if (flt) continue;
                                ++flt;
                                fprintf(stderr, "n0: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                if (!(ctgs->seq[ctg_idx].aux & 0x2)) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                            }  //precessor
                        } else { //potential to be a precessor 
                            /*if (get_max(c->cnts[j].cnt, 4) < 2) {*/
                            /*if (get_fmax(density, 4) < 2) {*/
                            /*if (tl < hd || get_max(c->cnts[j].cnt, 4) < 2) {*/
                            if (tl < hd || norm_cdf(tl, 0.5, tl + hd) <= 0.95) {
                                if (flt) continue;
                                ++flt;
                                if (ctg_idx2 + 1 != ctg_idx) {
                                    fprintf(stderr, "m1: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                    if (!(ctgs->seq[ctg_idx].aux & 0x2)) {
                                        if (norm_cdf(hd, 0.5, tl + hd) > .95) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                                        else susp_hd = j;
                                    }
                                } 
                            } // a predecessor
                            else {
                                if (fgt) continue;
                                ++fgt;
                                fprintf(stderr, "m0: %s\t%s\n", name, ctgs->seq[ctg_idx2].name);
                                kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                            } // a successor 
                        }
                        //validate 
                }
            }
            /*if (fgt && flt) break;*/
        }
		
        if (susp_hd != 0xFFFFFFFF) {
            int insert = 1;
            for (j = susp_hd + 1; j < susp_hd + 2 && j < c->n_cnt; ++j) {
                uint32_t sum = 0;
                for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
                uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
                uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (scf_idx2 == scf_idx && ctg_idx2 + 1 == ctg_idx) {
                fprintf(stderr, "val:%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                    if (hd > tl && norm_cdf(hd, 0.5, tl + hd) > 0.95) 
                        insert = 0;
                    break;   
                } 
            }
            if (insert) kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);      
        }
        if (susp_tl != 0xFFFFFFFF) {
            int insert = 1;
            for (j = susp_tl + 1; j < susp_tl + 2 && j < c->n_cnt; ++j) {
                uint32_t sum = 0;
                for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
                uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
                uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (scf_idx2 == scf_idx && ctg_idx2 == ctg_idx + 1) {
                fprintf(stderr, "val:%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                    if (hd < tl && norm_cdf(tl, 0.5, tl + hd) >= 0.95) 
                        insert = 0;
                    break;   
                } 
            }
            if (insert) kv_push(uint32_t, v, (ctg_idx) << 1 | 1);      
        }
        if (fgt == 0 && (ctgs->seq[ctg_idx].aux & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].aux >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        if (ctgs->seq[ctg_idx].aux & 1) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
	}
    *n_brks = v.n;
	return v.a;
}
uint32_t *find_breaks(cdict2_t *cds, sdict_t *ctgs, uint32_t *n_brks)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	for ( i = 0; i < n_cds; ++i) {
        uint32_t scf_idx = ctgs->seq[i].le;
        char *name = ctgs->seq[i].name;
        uint32_t ctg_idx = i; 
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j, z;
		c = cds + i;
        uint32_t fgt, flt;
        fgt = flt = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
            uint32_t sum = 0;
            for ( z = 0; z < 4; ++z) sum += c->cnts[j].cnt[z]; 
            uint32_t ctg_idx2 = sd_get(ctgs, c->cnts[j].name);
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
                //is_suc 
                fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);
                ctgn += 1;
                    if (ctgn >= 10) break;
                uint32_t hh, ht, th, tt, tl, hd, tl2, hd2;
                hh = c->cnts[j].cnt[0], ht = c->cnts[j].cnt[1], th = c->cnts[j].cnt[2], tt = c->cnts[j].cnt[3]; 
                tl = th + tt;
                hd = hh + ht;
                hd2 = hh + th;
                tl2 = ht + tt;
                if (ctg_idx2 > ctg_idx) {
                    /*fprintf(stderr, "suc: %s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                    if (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
						fprintf(stderr, "%u\t%u\n%s:n1\n", ctg_idx2, ctg_idx, name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } else { // validate direction
                        if (!(tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95))  {
                            /*fprintf(stderr, "%s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                            fprintf(stderr, "%s:m1\n", name);
                            kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                        }
                    }
                } else if (ctg_idx2 < ctg_idx) {
                    /*fprintf(stderr, "pre: %s\t%s\t%u\t%u\t%u\t%u\n", name, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        fprintf(stderr, "%s:n0\n", name);
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } else {
                        if (!((hd > tl || norm_cdf(tl, 0.5, tl + hd) <= 0.95) && (tl2 > hd2 || norm_cdf(hd2, 0.5, tl2 + hd2) <= 0.95))) {
                            fprintf(stderr, "%s:m0\n", name);
                            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                        }
                    }
                }
                /*
                if (ctg_idx2 > ctg_idx) {
                    if (tl < hd) continue; 
                    if (!fgt) ++fgt;
                    else continue;
                    if (ctg_idx2 != ctg_idx + 1) {
						fprintf(stderr, "%s:n1\n", name);
                        kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                    } else { // validate direction
                        if (!((tl > hd || norm_cdf(hd, 0.5, tl + hd) <= 0.95) && (hd2 > tl2 || norm_cdf(tl2, 0.5, tl2 + hd2) <= 0.95)))  {
							fprintf(stderr, "%s:m1\n", name);
                            kv_push(uint32_t, v, ctg_idx << 1 | 1);         
                        }
                    }
                } else {
                    if (tl > hd) continue; 
                    if (!flt) ++flt;
                    else  continue;
                    if (ctg_idx2 + 1 != ctg_idx) {
                        fprintf(stderr, "%s:n0\n", name);
                        kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                    } else {
                        uint32_t tl, hd, tl2, hd2;
                        tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];
                        hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];
                        hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];
                        tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];
                        if (!((hd > tl || norm_cdf(tl, 0.5, tl + hd) <= 0.95) && (tl2 > hd2 || norm_cdf(hd2, 0.5, tl2 + hd2) <= 0.95))) {
							fprintf(stderr, "%s:m0\n", name);
                            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);         
                        }
                    }
                } 
            */
            }
            if (fgt && flt) break;
		}
        if (fgt == 0 && (ctgs->seq[ctg_idx].aux & 1) == 0) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
        if (flt == 0 && (ctgs->seq[ctg_idx].aux >> 1 & 1) == 0)  
            kv_push(uint32_t, v, (ctg_idx - 1) << 1 | 1);  
        if (ctgs->seq[ctg_idx].aux & 1) 
            kv_push(uint32_t, v, ctg_idx << 1 | 1);  
	}
    *n_brks = v.n;
	return v.a;
}

int cmp_brks(const void *a, const void *b)
{
    uint32_t p = *(uint32_t *)a;
    uint32_t q = *(uint32_t *)b;
    if (p < q) return -1;
    else if (p == q) return 0;
    else return 1;
}

typedef struct {
	char *ctgn;
	uint32_t s, e, nl;
}name_t;

int parse_name(char *s, int l, name_t *nt)
{
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != ':' && s[i] != '-') continue;
		s[i] = 0;
		if (t == 0) nt->ctgn = q, nt->nl = i;
		else if (t == 1) nt->s = strtol(q, &r, 10), s[i] = '-';
		else if (t == 2) nt->e = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;

}

int sel_supt_reg(uint32_t *brks, uint32_t n_brks, sdict_t *ctgs) 
{
    qsort(brks, n_brks, sizeof(uint32_t), cmp_brks);  
    
    uint32_t i, j;
    /*for (i = 0; i < n_brks; ++i) {*/
        /*fprintf(stdout, "%s\t%d\t%d\n", ctgs->seq[brks[i]>>1].name, brks[i]>>1, brks[i]&1);*/
    /*}*/
    uint32_t stp = 0;
	fprintf(stdout, "track name=\"%s\" description=\"%s\"\n", "HiC", "HiC coverage");	
    name_t nt, nt1;
    for ( i = 1, j = 0; i <= n_brks; ++i) {
        if (i == n_brks || brks[i] != brks[j]) {
            //output  
            char *name = ctgs->seq[stp].name;
            parse_name(name, strlen(name), &nt); 
            uint32_t ctg_idx = brks[j];
            if (ctg_idx & 1) { // to the end
                /*fprintf(stdout, "%s\n", ctgs->seq[ctg_idx >> 1].name);*/
                if (ctg_idx >> 1 != stp) parse_name(ctgs->seq[ ctg_idx >> 1].name, strlen(ctgs->seq[ctg_idx>>1].name), &nt1); else nt1 = nt;
                fprintf(stdout, "%s\t%u\t%u\n", nt.ctgn, nt.s - 1, nt1.e);
                name[nt.nl] = ':';             
                ctgs->seq[ctg_idx >> 1].name[nt.nl] = ':';             
                stp = (ctg_idx>>1) + 1 ;
            } else {
                if ((ctg_idx >> 1) != stp + 1)parse_name(ctgs->seq[(ctg_idx >> 1)- 1].name, strlen(ctgs->seq[(ctg_idx >> 1) - 1].name), &nt1); else nt1 = nt;
                fprintf(stdout, "%s\t%u\t%u\n", nt.ctgn, nt.s - 1, nt1.e);
                name[nt.nl] = ':';             
                ctgs->seq[(ctg_idx>>1) - 1].name[nt.nl] = ':';             
                stp = (brks[i]>>1);
            } 
            j = i;        
        } 
    }
    return 0;
}

int ast_hic(char *fn, char *edge_fn, int min_wt, int use_sat, char *out_fn)
{
	sdict_t *ctgs = 0;
    if (!fn) {
        fprintf(stderr, "[E::%s] please set reference index file with -c\n", __func__);
        return 1;
    }
#ifdef VERBOSE
fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
    ctgs = col_ctgs(fn);	

	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict2_t* cds = calloc(n_ctg, sizeof(cdict2_t)); 
	uint32_t i;
	for ( i = 0; i < n_ctg; ++i) cd2_init(cds+i); 
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links(edge_fn, cds, ctgs);
    /*print_lnks(cds, ctgs);*/
	/*anothernorm(cds, ctgs);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] normalize contact matrix\n", __func__);
#endif
    norm_links(cds, ctgs);
	/*return 0;*/
     //sort cd2
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] sort contact matrix\n", __func__);
#endif
	for ( i = 0; i < n_ctg; ++i) cd2_sort(cds+i); 
    /*print_lnks(cds, ctgs);*/
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] find breaks\n", __func__);
#endif
    uint32_t n_brks;
    uint32_t *brks = find_breaks7(cds, ctgs, &n_brks);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] output support regions\n", __func__);
#endif
    
    sel_supt_reg(brks, n_brks, ctgs); 
	for (i = 0; i < n_ctg; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd2_destroy(cds +i);	

	} 
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
    free(brks);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	return 0;

}

int main(int argc, char *argv[])
{
	int c;
	uint32_t min_wt = 5; char *program = argv[0];
	char *ctg_fn = 0, *out_fn = 0;
	int use_sat = 0;
	while (~(c=getopt(argc, argv, "w:o:s:c:h"))) {
		switch (c) {
			case 'w': 
				min_wt = atoi(optarg);
				break;
			case 'c': 
				ctg_fn = optarg;
				break;
			case 'o': 
				out_fn = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s [<options>] <REF_FAI> <LINKS_MATRIX> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -w    INT      minimum weight for links [5]\n");
				fprintf(stderr, "         -s    FILE     sat file [nul]\n");
				fprintf(stderr, "         -o    FILE     output file [stdout]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require contact matrix and reference index file!\n", __func__); goto help;
	}
    char *sat_fn = argv[optind++];
	char *lnk_fn = argv[optind];
	fprintf(stderr, "Program starts\n");	
	int ret;
	ret = ast_hic(sat_fn, lnk_fn, min_wt, 0, out_fn);

	fprintf(stderr, "Program ends\n");	
	return ret;	

}

