CC      =  gcc
CFLAGS  =  -g -Wall -D VERBOSE -O2 -D PRINT_COVERAGE 
LDFLAGS = -lz -lm

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = ast_pb ast_10x ast_bion acc  pchlst detgaps ast_hic2 col_conts union_brks ast_bion_bnx split_fa # ast_postproc

.SUFFIXS:.c .o

all:$(PROG)
	mkdir -p ../bin && cp $(PROG) ../bin

ast_pb: ast_pb.o ast.o sdict.o paf.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

ast_10x: ast_10x.o ast.o sdict.o bamlite.o bed.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

col_conts: col_hic_lnks.o bamlite.o bed.o sdict.o cdict.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#hicbam2bed: hic_bam2bed.o bamlite.o bed.o sdict.o cdict.o 
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#ast_hic: ast_hic.o ast.o sdict.o bamlite.o  bed.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

ast_hic2: ast_hic2.o sdict.o  bed.o cdict.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

split_fa: split_fa.o sdict.o  bed.o cdict.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#ast_hic_ext: ast_hic_ext.o ast.o sdict.o bamlite.o  bed.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#ast_hic_ext2: ast_hic_ext2.o ast.o sdict.o bamlite.o  bed.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#ast_hic_diag: ast_hic_diag.o ast.o sdict.o bamlite.o  bed.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
ast_bion: ast_bion.o amap.o sdict.o ast.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#ast_bion_bnx_ref: ast_bion_bnx_ref.o amap.o sdict.o ast.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

ast_bion_bnx: ast_bion_bnx.o amap.o sdict.o ast.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
detgaps: det_gaps.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
#ast_val: ast_val.o bed.o sdict.o ast.o amap.o 
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


acc: acc.o sdict.o bed.o ast.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#breaks: breaks.o sdict.o bed.o ast.o
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


#join:joinh.o sdict.o bamlite.o 
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

#union:union_beds.o sdict.o bed.o ast.o 
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

union_brks:union_brks.o sdict.o bed.o ast.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
pchlst: gen_pchlst.o bed.o ast.o sdict.o bed.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
#ast_postproc: ast_postproc.o bed.o ast.o paf.o sdict.o 
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
#get_seqs: get_seqs.o sdict.o 
	#$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

