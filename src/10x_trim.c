/*
 * =====================================================================================
 *
 *       Filename:  10x.c
 *
 *    Description:  process 10x data to extract barcodes
 *
 *        Version:  1.0
 *        Created:  28/05/2018 14:38:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dg539@cam.ac.uk
 *   Organization:  Department of Genetics, Cambridge University
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>
#include <getopt.h>
#include "bseq.h"


#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef struct {
	int b_sz:30, fmt:1, cmp:1; // batch size and output format
	int bc_len, trimmed;
	char *out_prfx; //output prefix
	char *out_dir;	
}opts;

typedef struct {
	char *s;
	int	  l,m;
}buf_t;
/**
 * 
 * @func 
 * @
 *
 *
 */
int out_fa(eg_bseq1_t *s, eg_bseq1_t *r, FILE  *fp1, FILE *fp2, int bc_len, int trimmed)
{
	FILE *t = fp1;
	int jump = trimmed;
	eg_bseq1_t *rs = s;

	int i;
		
	for ( i = 0; i < 2; ++i) {
		fputc('>', t); fputs(rs->name, t); fputc('_', t);rs->seq[bc_len] = 0; fputs(rs->seq, t); fputc('\n', t); //first line
		fputs(rs->seq + jump, t); //sec line
		fputc('\n',t);	
		t = fp2, jump = 0, rs = r;
	} 
	return 0;
}
int out_fa2gz(eg_bseq1_t *s, eg_bseq1_t *r, gzFile  *fp1, gzFile *fp2, buf_t *b, int bc_len, int trimmed)
{
	gzFile *t = fp1;
	int jump = trimmed;
	eg_bseq1_t *rs = s;
	
	int l_name = strlen(rs->name);
	int buf_len = 4 + bc_len + l_name + (r->l_seq > s->l_seq ? r->l_seq : s->l_seq) + 32;	
	if (!(b->s && b->m >= buf_len)) {
		b->m = buf_len;
		b->s = realloc(b->s, sizeof(char) * kroundup32(b->m));
	}
	
	int i;
		
	for ( i = 0; i < 2; ++i) {
		int j = 0;
		b->s[j++] = '>', strcpy(b->s + j,s->name), j += l_name;
		b->s[j++] = '_', strncpy(b->s + j, s->seq, bc_len), j += bc_len; //first line
		b->s[j++] = '\n', strcpy(b->s + j, rs->seq+jump), j += rs->l_seq - jump;
		b->s[j++] = '\n', b->s[j] = 0;
		gzwrite(t, b->s, j);
		t = fp2, jump = 0, rs = r;
	} 
	return 0;
}

int out_fq2gz(eg_bseq1_t *s, eg_bseq1_t *r, gzFile *fp1, gzFile *fp2, buf_t *b, int bc_len, int trimmed)
{

	int l_name = strlen(s->name);
	int buf_len = 7 + bc_len + l_name + ((r->l_seq > s->l_seq ? r->l_seq : s->l_seq) << 1) + 32;	
	if (!(b->s && b->m >= buf_len)) {
		b->m = buf_len;
		b->s = realloc(b->s, sizeof(char) * kroundup32(b->m));
	}
	
	gzFile *t = fp1;
	eg_bseq1_t *rs = s;
	int jump = trimmed;
	int i;	
	for ( i = 0; i < 2; ++i) {
		int j = 0;
		b->s[j++] = '@', strcpy(b->s + j,s->name), j += l_name;
		b->s[j++] = '_', strncpy(b->s + j, s->seq, bc_len), j += bc_len; //name should be the same
		b->s[j++] = '\n', strcpy(b->s + j, rs->seq+jump), j += rs->l_seq - jump;
		strcpy(b->s+j, "\n+\n"), j += 3;
		strcpy(b->s+j, rs->qual + jump), j += rs->l_seq - jump;
		b->s[j++]='\n', b->s[j] = 0;
		gzwrite(t, b->s, j);
		t = fp2, rs = r,jump = 0;	
	} 
	return 0;
}

int out_fq(eg_bseq1_t *s, eg_bseq1_t *r, FILE *fp1, FILE *fp2, int bc_len, int trimmed)
{
	FILE *t = fp1;
	eg_bseq1_t *rs = s;
	int jump = trimmed;
	int i;	
	for ( i = 0; i < 2; ++i) {
		fputc('@', t); fputs(s->name, t); fputc('_', t);s->seq[bc_len] = 0; fputs(s->seq, t); fputc('\n', t); //first line
		fputs(rs->seq + jump, t); //sec line
		fputs("\n+\n", t); // third
		fputs(rs->qual + jump, t); //forth
		fputc('\n',t);	
		t = fp2, jump = 0, rs = r;	
	} 
	return 0;
}

int dump_reads2gz(eg_bseq1_t *egs, int n_seqs, gzFile *fp1, gzFile *fp2, int fmt, int bc_len, int trimmed)
{
	int i;
	buf_t *b = calloc(1, sizeof(buf_t));
	for (i=0; i < n_seqs; i += 2) {
		//output sequence name
		if (fmt) out_fq2gz(egs + i, egs+i+1, fp1, fp2, b, bc_len, trimmed);
		else out_fa2gz(egs+i, egs+i+1, fp1, fp2, b, bc_len, trimmed);	
	}
	if (b && b->s) {free(b->s); free(b);}
	return 0;
}

int dump_reads(eg_bseq1_t *egs, int n_seqs, FILE *fp1, FILE *fp2, int fmt, int bc_len, int trimmed)
{
	int i;
	for (i=0; i < n_seqs; i += 2) {
		//output sequence name
		if (fmt) out_fq(egs + i, egs+i+1, fp1, fp2, bc_len, trimmed);
		else out_fa(egs+i, egs+i+1, fp1, fp2, bc_len, trimmed);	
	}
	return 0;
}

int proc_10x_readsgz(char *fn[], int n_fn, opts *o)
{
		
	eg_bseq_file_t **fps = (eg_bseq_file_t **)(calloc(n_fn, sizeof(eg_bseq_file_t *))); 
	
	int i;
	if (n_fn != 2) {
		fprintf(stderr, "[E::%s] not pairwise files\n", __func__);
	   	return 1; 
	}
	for ( i = 0; i < n_fn; ++i) {
		fps[i] = eg_bseq_open(fn[i]);
		if (fps[i] == 0) {
			fprintf(stderr, "[E::%s] fail to open %s\n", __func__, fn[i]);
			int j; 
			for ( j = 0; j < i; ++j ) eg_bseq_close(fps[j]);
			if (fps) free(fps);
			return 1;
		}
	}
	//open output file handles
	int  l_dir = 1; char *dir = "."; 
	if (o->out_dir) l_dir = strlen(o->out_dir), dir = o->out_dir;
	
	int l_prfx = 3; char *prefx = "out";
	if (o->out_prfx) l_prfx = strlen(o->out_prfx), prefx = o->out_prfx;
	char *fn_t = (char *)malloc(sizeof(char)*(l_dir + l_prfx + 12)); 
	sprintf(fn_t, "%s/%s_%c.%s", dir, prefx, '1', o->fmt ? "fq.gz" : "fa.gz");	
	fn_t[l_prfx + l_dir + 9] = 0;
	
	gzFile *fp1, *fp2;
	fp1 = gzopen(fn_t, "wb");
	fn_t[l_dir + l_prfx+2] = '2';
	fp2 = gzopen(fn_t, "wb");
	free(fn_t);
	if (!(fp1 && fp2)) {
		fprintf(stderr, "[E::%s] fail to open output files", __func__);
		return 1;
	}
	//read sequences
	eg_bseq1_t *seqs;
	while (1) {
		int n_seqs;
		seqs = eg_bseq_read_frag(n_fn, fps, o->b_sz, 1, &n_seqs);
		if (seqs) {
			dump_reads2gz(seqs, n_seqs, fp1, fp2, !!o->fmt, o->bc_len, o->trimmed); 	
			eg_bseq_destroy(seqs, n_seqs);	
		} else 
			break; //finished
	}
	gzclose(fp1);
	gzclose(fp2);	
	if (fps) {
		for ( i = 0; i < n_fn; ++i) eg_bseq_close(fps[i]);
		free(fps);
	}
	return 0;	
}


int proc_10x_reads(char *fn[], int n_fn, opts *o)
{
		
	eg_bseq_file_t **fps = (eg_bseq_file_t **)(calloc(n_fn, sizeof(eg_bseq_file_t *))); 
	
	int i;
	if (n_fn != 2) {
		fprintf(stderr, "[E::%s] not pairwise files\n", __func__);
	   	return 1; 
	}
	for ( i = 0; i < n_fn; ++i) {
		fps[i] = eg_bseq_open(fn[i]);
		if (fps[i] == 0) {
			fprintf(stderr, "[E::%s] fail to open %s\n", __func__, fn[i]);
			int j; 
			for ( j = 0; j < i; ++j ) eg_bseq_close(fps[j]);
			if (fps) free(fps);
			return 1;
		}
	}
	//open output file handles
	int  l_dir = 1; char *dir = "."; 
	if (o->out_dir) l_dir = strlen(o->out_dir), dir = o->out_dir;

	int l_prfx = 3; char *prefx = "out";
	if (o->out_prfx) l_prfx = strlen(o->out_prfx), prefx = o->out_prfx;

	char *fn_t = (char *)malloc(sizeof(char)*(l_dir+l_prfx + 8)); 
	sprintf(fn_t, "%s/%s_%c.%s", dir, prefx, '1', o->fmt ? "fq" : "fa");	
	fn_t[l_prfx + l_dir +  6] = 0;
	
	FILE *fp1, *fp2;
	fp1 = fopen(fn_t, "w");
	fn_t[l_dir + l_prfx+2] = '2';
	fp2 = fopen(fn_t, "w");
	free(fn_t);
	if (!(fp1 && fp2)) {
		fprintf(stderr, "[E::%s] fail to open output files", __func__);
		return 1;
	}
	//read sequences
	eg_bseq1_t *seqs;
	while (1) {
		int n_seqs;
		seqs = eg_bseq_read_frag(n_fn, fps, o->b_sz, 1, &n_seqs);
		if (seqs) {
			dump_reads(seqs, n_seqs, fp1, fp2, !!o->fmt, o->bc_len, o->trimmed); 	
			eg_bseq_destroy(seqs, n_seqs);	
		} else 
			break; //finished
	}
	fclose(fp1);
	fclose(fp2);	
	if (fps) {
		for ( i = 0; i < n_fn; ++i) eg_bseq_close(fps[i]);
		free(fps);
	}
	return 0;	
}





int main(int argc, char *argv[]) 
{
	int c;
	opts o;
	o.fmt = 1;
	o.b_sz = 100000000;//100M default
	o.out_prfx = NULL;
	o.cmp = 0;
	o.trimmed = 23;
	o.bc_len = 16;
	o.out_dir = NULL; 
	while (~(c=getopt(argc, argv, "F:b:o:p:l:m:ch"))) {
		switch(c) {
			case 'F':
				if (!strncmp(optarg, "FASTQ", 5))
					o.fmt = 1;
				else 
					o.fmt = 0;
				break;
			case 'b':
				o.b_sz = atoi(optarg);	// only 31 bits be careful
				break;	
			case 'l':
				o.bc_len = atoi(optarg);	// only 31 bits be careful
				break;	
			case 'm':
				o.trimmed = atoi(optarg);	// only 31 bits be careful
				break;	
			case 'p':
				o.out_prfx = optarg;
				break;
			case 'o':
				o.out_dir = optarg;
				break;
			case 'c':
				o.cmp = 1;
				break;
			default:
				if (c != 'h') fprintf(stderr, "\n[E::%s] undefined option %c\n", __func__, c);
help:
				fprintf(stderr, "\nUsage: 10x [options...] <in1.fq/in1.fq.gz> <in2.fq/in2.fq.gz>\n\n");
				fprintf(stderr, "Options:   -b    INT    batch size [100M]\n");
				fprintf(stderr, "           -l    INT    barcode length [16]\n");
				fprintf(stderr, "           -m    INT    trimmed off bases [23]\n");
				fprintf(stderr, "           -p    STR    prefix of output files [out]\n");
				fprintf(stderr, "           -o    STR    output direcotory [.]\n");
				fprintf(stderr, "           -F    STR    output file format: FASTA or FASTQ [FASTQ]\n");
				fprintf(stderr, "           -c    BOOL   compressed [false]\n");
				fprintf(stderr, "           -h           help\n");
				fprintf(stderr, "\n");
				return 1;	
		}
	}
	int n_fn = argc - optind;
		
	if (n_fn) {
		if (o.cmp) {
			if (!proc_10x_readsgz(argv+optind, n_fn, &o)) fprintf(stderr, "finish process 10x sequences\n"); else goto help;
		} else {
			if (!proc_10x_reads(argv+optind, n_fn, &o)) fprintf(stderr, "finish process 10x sequences\n"); else goto help;
		}
	} else {
		fprintf(stderr, "[E::%s] original pairwise fastq files required!\n", __func__);
		goto help;
	}
}




