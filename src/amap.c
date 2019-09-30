/*
 * =====================================================================================
 *
 *       Filename:  amap.c
 *
 *    Description:  read bionano map file
 *
 *        Version:  1.0
 *        Created:  17/09/2018 09:24:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include "amap.h"

#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)

amap_file_t *amap_open(const char *fn)
{
	kstream_t *ks;
	gzFile fp;
	amap_file_t *amp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	amp = (amap_file_t*)calloc(1, sizeof(amap_file_t));
	amp->fp = ks;
	return amp;
}

int amap_close(amap_file_t *amp)
{
	kstream_t *ks;
	if (amp == 0) return 0;
	free(amp->buf.s);
	ks = (kstream_t*)amp->fp;
	gzclose(ks->f);
	ks_destroy(ks);
	free(amp);
	return 0;
}

int amp_parse_rqmap(int l, char *s, rqmap_entry *pr) // s must be NULL terminated
{ // on return: <0 for failure; 0 for success; >0 for filtered
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->cmap_id = strtol(q, &r, 10);
		else if (t == 1) pr->ctg_len = strtol(q, &r, 10);
		else if (t == 2) pr->num_sites = strtol(q, &r, 10);
		else if (t == 3) pr->site_id = strtol(q, &r, 10);
		else if (t == 4) pr->label_chan = strtol(q, &r, 10);
		else if (t == 5) pr->pos = strtof(q, &r); 
		else if (t == 6) pr->std_dev = strtof(q, &r);
		else if (t == 7) pr->cov = strtof(q, &r);
		/*else if (t == 8) pr->te = strtol(q, &r, 10);*/
		/*else if (t == 9) pr->ml = strtol(q, &r, 10);*/
		/*else if (t == 10) pr->bl = strtol(q, &r, 10);*/
		/*else if (t == 11) pr->mq = strtol(q, &r, 10);*/
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 7) return -1;
	return 0;
}
int amp_parse_k(int l, char *s, key_entry *pr) // s must be NULL terminated
{ // on return: <0 for failure; 0 for success; >0 for filtered
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->ctg_id = strtol(q, &r, 10);
		else if (t == 1) {
			//add code to trim off chars after space
			char* w;
			for (w = q; w < s + i; ++w) if (*w == ' ') {*w = 0; break;} 
			pr->ctg_name = q;
		} 
		else if (t == 2) pr->ctg_len = strtof(q, &r);
		/*else if (t == 3) pr->qry_s = strtof(q, &r);*/
		/*else if (t == 4) pr->qry_e = strtof(q, &r);*/
		/*else if (t == 5) pr->ref_s = strtof(q, &r); */
		/*else if (t == 6) pr->ref_e = strtof(q, &r);*/
		/*else if (t == 7) pr->ori = q[0];*/
		/*else if (t == 8) pr->conf = strtof(q, &r);*/
		/*else if (t == 9) pr->hit_enum = q;*/
		/*else if (t == 10) pr->qry_len = strtof(q, &r);*/
		/*else if (t == 11) pr->ref_len = strtof(q, &r);*/
		/*else if (t == 12) pr->label_chan = strtol(q, &r, 10);*/
		/*else if (t == 13) pr->alignment = q;*/
		/*else if (t == 8) pr->te = strtol(q, &r, 10);*/
		/*else if (t == 9) pr->ml = strtol(q, &r, 10);*/
		/*else if (t == 10) pr->bl = strtol(q, &r, 10);*/
		/*else if (t == 11) pr->mq = strtol(q, &r, 10);*/
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;
}

int amp_parse_x(int l, char *s, xmap_entry *pr) // s must be NULL terminated
{ // on return: <0 for failure; 0 for success; >0 for filtered
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->entry_id = strtol(q, &r, 10);
		else if (t == 1) pr->qry_id = strtol(q, &r, 10);
		else if (t == 2) pr->ref_id = strtol(q, &r, 10);
		else if (t == 3) pr->qry_s = strtof(q, &r);
		else if (t == 4) pr->qry_e = strtof(q, &r);
		else if (t == 5) pr->ref_s = strtof(q, &r); 
		else if (t == 6) pr->ref_e = strtof(q, &r);
		else if (t == 7) pr->ori = q[0];
		else if (t == 8) pr->conf = strtof(q, &r);
		else if (t == 9) pr->hit_enum = q;
		else if (t == 10) pr->qry_len = strtof(q, &r);
		else if (t == 11) pr->ref_len = strtof(q, &r);
		else if (t == 12) pr->label_chan = strtol(q, &r, 10);
		else if (t == 13) pr->alignment = q;
		/*else if (t == 8) pr->te = strtol(q, &r, 10);*/
		/*else if (t == 9) pr->ml = strtol(q, &r, 10);*/
		/*else if (t == 10) pr->bl = strtol(q, &r, 10);*/
		/*else if (t == 11) pr->mq = strtol(q, &r, 10);*/
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 13) return -1;
	return 0;
}
int amap_read_rq(amap_file_t *pf, rqmap_entry *r)
{
	int ret, dret;
file_read_more:
	ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
	/*fprintf(stderr, "%s", pf->buf.s);*/
	if (ret < 0) return ret;
	if (pf->buf.s[0] == '#') 
		goto file_read_more; //skip comments 
	else 
		ret = amp_parse_rqmap(pf->buf.l, pf->buf.s, r);
	if (ret < 0) goto file_read_more;
	return ret;
}
int amap_read_x(amap_file_t *pf, xmap_entry *r)
{
	int ret, dret;
file_read_more:
	ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
	if (ret < 0) return ret;
	if (pf->buf.s[0] == '#') 
		goto file_read_more; //skip comments 
	else 
		ret = amp_parse_x(pf->buf.l, pf->buf.s, r);
	if (ret < 0) goto file_read_more;
	return ret;
}
int amap_read_k(amap_file_t *pf, key_entry *r)
{
	int ret, dret;
file_read_more:
	ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
	
	if (ret < 0) return ret; //should be finished
	if (pf->buf.s[0] == '#') 
		goto file_read_more; //skip comments 
	else 
		ret = amp_parse_k(pf->buf.l, pf->buf.s, r);
	if (ret < 0) goto file_read_more;
	return ret;
}
