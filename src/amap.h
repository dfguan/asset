/*
 * =====================================================================================
 *
 *       Filename:  amap.h
 *
 *    Description:  for reference query cmap and xmap 
 *
 *        Version:  1.0
 *        Created:  16/09/2018 20:09:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef AMAP_H
#define AMAP_H

#include <stdint.h>
#include <sys/types.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

typedef struct {
	void *fp;
	kstring_t buf;
} amap_file_t;

typedef struct {
	const char *qn, *tn; // these point to the input string; NOT allocated
	uint32_t ql, qs, qe, tl, ts, te;
	uint32_t ml:31, rev:1, bl;
	int mq;
} paf_rec_t;

typedef struct {
	int cmap_id;
	float ctg_len;
	int		num_sites;
	int		site_id;
	int		label_chan;
	float	pos, std_dev, cov, occ;
}rqmap_entry;

typedef struct {
	int	entry_id, qry_id, ref_id;
	float qry_s, qry_e, ref_s, ref_e;	
	float conf;
	char ori;
	char* hit_enum;
	float qry_len, ref_len;
	int  label_chan;
	char *alignment;
}xmap_entry;

typedef struct {
	int		ctg_id;
	char	*ctg_name;
	float	ctg_len;
}key_entry;


#ifdef __cplusplus
extern "C" {
#endif

amap_file_t *amap_open(const char *fn);
int amap_close(amap_file_t *pf);
int amap_read_rq(amap_file_t *pf, rqmap_entry *r);
int amap_read_x(amap_file_t *pf, xmap_entry *r);
int amap_read_k(amap_file_t *pf, key_entry *r);

#ifdef __cplusplus
}
#endif

#endif
