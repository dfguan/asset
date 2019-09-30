/*
 * =====================================================================================
 *
 *       Filename:  det_gaps.c
 *
 *    Description:  detect gaps
 *
 *        Version:  1.0
 *        Created:  14/10/2018 11:42:57
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
#include <stdint.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread, gzseek)

int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;
  int l;
  char *program; 	
  (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.fasta>\n", program);
    return 1;
  }
  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
	uint32_t length = seq->seq.l;
	uint32_t i;
	if (!(seq->seq.s[0] == 'n' || seq->seq.s[0] == 'N'))
		fprintf(stdout, "%s\t%u\t%u\n", seq->name.s, 0, 0);
	for ( i = 0; i < length; ++i) {
		if (seq->seq.s[i] == 'n' || seq->seq.s[i] =='N') {
			uint32_t j;
			for ( j = i + 1; j <= length; ++j) {
				if (j == length ||(seq->seq.s[j] != 'n' && seq->seq.s[j] != 'N')) { 
					fprintf(stdout, "%s\t%u\t%u\n", seq->name.s, i, j);
					break;
				}	
			}
			i = j;
		}		
	}
	fprintf(stdout, "%s\t%u\t%u\n", seq->name.s, length, length);
  }
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}


