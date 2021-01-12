/*
 * =====================================================================================
 *
 *       Filename:  aa.h
 *
 *    Description:  header for aa file
 *
 *        Version:  1.0
 *        Created:  15/09/2018 19:13:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef AA_H
#define AA_H

#include <stdio.h>
#include <stdint.h>
#include "sdict.h"

#define CONT_THRES 20
#define CONT_THRES_10X 5

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
	//int	*p;//need to change for future use?
	uint32_t *p;
	int n,m;	
}pos_t;

typedef struct {
	pos_t	*ctg_pos;
	int		n,m; //contig number
}ctg_pos_t;

typedef struct {
	uint32_t s, e;
	int coverage;
}cov_t;

typedef struct {
	uint32_t s, e;
}cors;

typedef struct {
	cors * coords;
	int sorted;
	char **anno;
	size_t n, m;
}cord_t;

typedef struct {
	int		n,m;
	uint64_t len, tot_cov;
	cov_t *intv;
}cov_ary_t;


typedef struct {
	cord_t *ct;
	size_t n, m;	
}ns_t;

#ifdef __cplusplus
extern "C" {
#endif

void cov_ary_push(cov_ary_t *c, uint32_t s, uint32_t e, int cov);
void cov_ary_destroy(cov_ary_t *ca, int n);

ctg_pos_t *ctg_pos_init();
void ctg_pos_push(ctg_pos_t *_d,int s);
void ctg_pos_reset(ctg_pos_t *_d);
void ctg_pos_destroy(ctg_pos_t *_d);

void pos_push(pos_t *ps, uint32_t _p);
void pos_destroy(pos_t *ps);


void cord_push(cord_t *c, cors *cord);
void cord_push1(cord_t *c, cors *cord, char *ann);
int cord_sort(cord_t *c, int n);
int cord_merge(cord_t *c, int n);
void cord_destroy(cord_t *c, int n);
void print_coverage(cov_ary_t *ca, sdict_t* ctgs, char *tp);
void print_coverage_wig(cov_ary_t *ca, sdict_t* ctgs, char *tp, uint32_t ws, char *out_dir);
void sel_sup_reg(cov_ary_t *ca, int min_cov, int max_cov, sdict_t* ctgs, char *tp, char *desc);
cov_ary_t *cal_cov(ctg_pos_t *d, sdict_t* ctgs, float *avgc);


void sel_sup_reg_dyn(cov_ary_t *ca, float min_cov_rat, int min_cov, int max_cov, sdict_t* ctgs, char *tp, char *desc);
void print_base_coverage(cov_ary_t *ca, sdict_t* ctgs, char *tp, char *out_dir);
void print_coverage_stat(cov_ary_t *ca, sdict_t* ctgs, char *tp, char *out_dir);
void print_coverage_by_wind(cov_ary_t *ca, sdict_t* ctgs, char *tp, uint32_t ws, char *out_dir);
void ns_push(ns_t *ns, uint32_t n);
void ns_destroy(ns_t *ns);
uint32_t cal_n50(uint32_t *v, uint32_t n);
#ifdef __cplusplus
}
#endif
#endif
