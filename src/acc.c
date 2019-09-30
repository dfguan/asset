/*
 * =====================================================================================
 *
 *       Filename:  acc.c
 *
 *    Description:  accumulate all the envidence
 *
 *        Version:  1.0
 *        Created:  18/09/2018 12:59:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
// use uint8_t:4l to keep number of support technology and the highest bit used for a gap, 4h for technology track

#include <stdio.h>
#include <zlib.h>

#include "ast.h"
#include "bed.h"
#include "sdict.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)

char *tp[] = {"PB", "TX", "HC", "BN", "UN"}; //E2 E3 E4 for evidence or more 

int acc_evd(char *evd, uint8_t *sig, sdict_t *ctgs, uint64_t blk_l)
{
	gzFile fp = gzopen(evd, "r");
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open %s\n",__func__, evd);
		return -1;
	}
	kstream_t *ks = ks_init(fp);
	kstring_t buf = (kstring_t){0,0,0};  //don't forget to initiate this variable 
	int dret;
	uint32_t max_len = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		char *r, *ctg_n, *q;
		uint32_t s,e;
		int t, i;
		for (i = t = 0, q= buf.s; i <= buf.l;++i) {
			if (i < buf.l && buf.s[i] != '\t') continue;	
			buf.s[i] = 0; 
			if (t == 0) ctg_n = q;
			else if (t == 1) s = strtol(q, &r, 10);
			else if (t == 2) e = strtol(q, &r, 10);
			++t, q = i < buf.l ? &buf.s[i+1] : 0;	
			/*fprintf(stderr, "%s\n",q);*/
		}
		//initate signals 
		int ind = sd_get(ctgs, ctg_n);	
		uint32_t j;
		uint8_t *ps = sig + blk_l * ind;	
		/*fprintf(stderr, "%s\t%u\t%u\n", ctg_n, s, e);*/
		for (j = s-1; j <= e-1; ++j) {
			if (j & 1) ps[j>>1] += 1;
			else ps[j>>1] += 0x10;
		} 
		/*fprintf(stderr, "%s\t%u\t%u\n", ctg_n, s, e);*/
	}	
	ks_destroy(ks);
	gzclose(fp);
	return max_len;
}

int acc_evd2(char *evd, uint8_t *sig, sdict_t *ctgs, uint64_t* sig_ind)
{
	bed_hdr_t hdr = (bed_hdr_t){"UN","unkown data type"};
	bed_file_t *bfp = bed_open(evd);
	if (!bfp) return -1;
	bed_hdr_read(bfp, &hdr);
	bed_rec_t r;
	uint8_t tp_s = sizeof(tp) / sizeof(char *);
	/*fprintf(stderr, "%hhu\n", tp_s);*/
	uint8_t ind_tp;
	for ( ind_tp = 0; ind_tp < tp_s; ++ind_tp) {
		if (*(uint16_t *)hdr.type == *(uint16_t *)tp[ind_tp]) break;
	}
		/*fprintf(stderr, "haha:%s\n", hdr.type);*/
	while (bed_read(bfp, &r) >= 0) {
		uint64_t ind = sd_get(ctgs, r.ctgn);
		
		uint8_t *ps = sig + sig_ind[ind];	
		uint32_t j;	
		for (j = r.s; j < r.e; ++j) {
			/*ps[];*/
			/*if (j & 1) ps[j>>1] += 1;*/
			/*else ps[j>>1] += 0x10;*/
			ps[j] += 1;
			ps[j] |= (ind_tp << 4); //only works when one evidence is there.
		} 
	}	
	bed_close(bfp);
	return 0;
}


uint64_t *init_sig_index(sdict_t *ctgs)
{
	uint64_t *a = NULL;
	if (ctgs->n_seq) {
		a = malloc(sizeof(uint64_t) * ctgs->n_seq);
		a[0] = 0;
		uint32_t i;
		for (i = 1; i < ctgs->n_seq; ++i) a[i] = a[i-1] + ctgs->seq[i-1].len;
	}
	return a;	
}


uint64_t col_ctgs(char *ctg_fn, sdict_t *ctgs)
{
	gzFile fp = gzopen(ctg_fn, "r");
	if (!fp) {
		fprintf(stderr, "[E::%s] fail to open %s\n",__func__, ctg_fn);
		return -1;
	}
	kstream_t *ks = ks_init(fp);
	kstring_t buf = (kstring_t){0,0,0} ;  
	int dret;
	uint64_t total_len = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		char *ctg_n = buf.s, *p, *r;
		for (p = buf.s; p < buf.s+buf.l && *p!='\t';++p);
		*p = 0;	
		uint32_t ctg_l = strtol(p+1, &r, 10);
		/*fprintf(stderr, "%s\t%u\n",ctg_n, ctg_l);*/
		sd_put(ctgs, ctg_n, ctg_l);
		/*if (max_len < ctg_l) */
		total_len += ctg_l;
	}	
	ks_destroy(ks);
	gzclose(fp);
	return total_len;
}

int gen_wig(uint8_t *signals, sdict_t *ctgs, uint32_t blk_len)
{
	uint32_t i;
	for (i = 0; i < ctgs->n_seq; ++i) { //for each contig
		uint8_t *p = signals + i * blk_len;
		uint32_t j = 0;
		uint32_t l_seq = ctgs->seq[i].len;
		uint32_t s[] = {0,0,0,0,0};
		while (j < l_seq) {
			uint32_t z;
			uint8_t d, pd = j & 1 ? p[j>>1] & 0x0F : p[j>>1] >> 4;
			for (z = j; z <= l_seq; ++z) {
				if (z < l_seq) d = z & 1 ? p[z>>1] & 0x0F : p[z>>1] >> 4;
				if (z ==l_seq || d != pd) {
					/*if (pd) {*/
						fprintf(stdout, "variableStep chrom=%s span=%u\n",ctgs->seq[i].name, z -j);	
						fprintf(stdout, "%u %hhu\n",j+1,pd);	
						s[pd] += z -j;
					/*}*/
					break;
				} 		
			}
			j = z;
			/*if (z == l_seq && d) {*/
				/*fprintf(stdout, "variableStep chrom=%s span=%u\n",ctgs->seq[i].name, z -j);	*/
				/*fprintf(stdout, "%u %hhu\n",j,d);	*/
			/*} //could be put in the loop? */
			/*j = z;*/
		}
		fprintf(stderr, "%s\t%u\t%u\t%u\t%u\n", ctgs->seq[i].name, s[1],s[2],s[3],s[4]);
	}	
	return 0;
}	

int gen_bed(uint8_t *signals, sdict_t *ctgs, uint64_t *sig_ind) 
{
	uint32_t i;
	for (i = 0; i < ctgs->n_seq; ++i) { //for each contig
		/*fprintf(stderr, "%d %s, %lu %lu\n", ctgs->n_seq, ctgs->seq[i].name, sig_ind[i], ctgs->seq[i].len);*/
		uint8_t *p = signals + sig_ind[i];
		uint32_t j = 0;
		uint32_t l_seq = ctgs->seq[i].len;
		uint32_t s[] = {0,0,0,0,0};
		uint32_t cnt[] = {0,0,0,0,0};
		while (j < l_seq) {
			uint32_t z;
			uint8_t ptp = (p[j] >> 4 ) & 0x0F, pd = p[j] & 0x0F;
			uint8_t tpi, d;
			for (z = j; z <= l_seq; ++z) {
				if (z < l_seq) 
					tpi = (p[z] >> 4 ) & 0x0F, d = p[z] & 0x0F;
				if (z ==l_seq || d != pd) {
					/*if (pd) {*/
						/*fprintf(stdout, "variableStep chrom=%s span=%u\n",ctgs->seq[i].name, z -j);	*/
						/*fprintf(stdout, "%u %hhu\n",j+1,pd);	*/
						if ((pd & 0x7) == 1) fprintf(stdout, "%s\t%u\t%u\t%hhu\t%s\t%c\n", ctgs->seq[i].name, j, z, pd & 0x7, tp[ptp], pd & 0x8 ? 'G':'N');
						else fprintf(stdout, "%s\t%u\t%u\t%hhu\t\t%c\n", ctgs->seq[i].name, j, z, pd & 0x7, pd & 0x8 ? 'G':'N');
					s[pd & 0x7] += z -j;
					cnt[pd & 0x7] += 1;
					/*}*/
					break;
				} else if ((pd & 0x7)== 1 && ptp != tpi) {
						fprintf(stdout, "%s\t%u\t%u\t%hhu\t%s\t%c\n", ctgs->seq[i].name, j, z, pd & 0x7, tp[ptp], pd & 0x8 ? 'G':'N');
						s[pd & 0x7] += z -j;
						cnt[pd & 0x7] += 1;
					/*}*/
					break;
				} 	
			}
			j = z;
			/*if (z == l_seq && d) {*/
				/*fprintf(stdout, "variableStep chrom=%s span=%u\n",ctgs->seq[i].name, z -j);	*/
				/*fprintf(stdout, "%u %hhu\n",j,d);	*/
			/*} //could be put in the loop? */
			/*j = z;*/
		}
		int q;
		fprintf(stderr, "%s", ctgs->seq[i].name);
		for (q = 0; q < 5; ++q){
			fprintf(stderr, "\t%lf:%u", (double)s[q]/ctgs->seq[i].len, cnt[q]);	
		}
		fprintf(stderr, "\n");	
	}	
	return 0;
}	

void init_gaps(char *gap_fn, ns_t *ns, sdict_t *ctgs)
{
	bed_file_t* bf = bed_open(gap_fn);
	bed_rec_t r;
	while (bed_read(bf, &r) >= 0) {
		uint32_t ind = sd_put(ctgs, r.ctgn, r.e);
		ctgs->seq[ind].len = r.e;
		/*if (r.e - r.s >= max_ins_len) {*/
		ns_push(ns, ind);
		cors tmp = (cors){r.s, r.e};
		cord_push(&ns->ct[ind], &tmp);				
		/*}*/
	}
	bed_close(bf);
}

uint64_t col_gs(sdict_t *sd)
{
	uint64_t tl = 0;
	int i;
	for ( i = 0; i < sd->n_seq; ++i) tl += sd->seq[i].len;
	return tl;
}

int mark_gaps(uint8_t* sigs, uint64_t *sig_idx, ns_t *ns)
{
	int i;	
	/*fprintf(stderr, "%d\n", ns->n);*/
	for ( i = 0; i < ns->n; ++i) {
		uint8_t *p = sigs + sig_idx[i];
		cord_t *ct = &ns->ct[i];	
		cors *cs = ct->coords;
		uint32_t j,z;
		for ( j = 0; j < ct->n; ++j) {
			/*fprintf(stderr, "%u\t%u\n", cs[j].s, cs[j].e);*/
			for ( z = cs[j].s; z < cs[j].e; ++z) 
				p[z] |= 0x8;		
		} 
	}	
	return 0;
}


int main(int argc, char *argv[])
{
	int n = argc;
		
	
	if (n < 2) {
		fprintf(stderr, "[E::%s] require at least one bed file\n", __func__);
		fprintf(stderr, "acc <GAP.BED> <BED_FILEs> ...\n");
		return -1;
	} 
	//initiate index based on ctgs name and length file 
										
	char *gap_fn = argv[1];
	ns_t *ns = calloc(1, sizeof(ns_t));
	sdict_t *ctgs = sd_init();
	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collect contig names and contig length\n", __func__);
#endif
	init_gaps(gap_fn, ns, ctgs);
	uint64_t total_len = col_gs(ctgs);
	/*fprintf(stderr, "%lu\n", total_len)	;*/
	/*uint32_t blk_len = (max_len & 1) + (max_len >> 1);	*/
	uint64_t *sig_ind = init_sig_index(ctgs);
	
	uint8_t *signals = calloc(total_len + 32, sizeof(uint8_t));
	
	mark_gaps(signals, sig_ind, ns);	
	ns_destroy(ns);
	if (!signals) {
		sd_destroy(ctgs);
		fprintf(stderr, "[E::%s] fail to allocate memory for signals, required %lu bytes\n", __func__, total_len + 32);	
		return -1;	
	}	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] accumulate evidence\n", __func__);
#endif
	int i;
	for (i = 2; i < argc; ++i) 
		/*int acc_evd2(char *evd, uint8_t *sig, sdict_t *ctgs, uint32_t blk_l)*/
		acc_evd2(argv[i], signals, ctgs, sig_ind);	

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] generate wig\n", __func__);
#endif
	/*gen_wig(signals, ctgs, blk_len);	*/
	/*fprintf(stderr, "before: %p\n", ctgs);*/
	gen_bed(signals, ctgs, sig_ind);	
	/*fprintf(stderr, "after: %p\n", ctgs);*/
		
	/*fprintf(stderr, "free contgis");*/
	sd_destroy(ctgs);
	/*fprintf(stderr, "free signals");*/
	free(signals);
	/*fprintf(stderr, "free signals index");*/
	if (sig_ind) free(sig_ind);
	return 0;
}

