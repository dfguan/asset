from runner.manager import manager
from runner.hpc import hpc
from multiprocessing import Process, Pool
import sys, os, json

# utilities
def checkf(p):
    return os.path.isfile(p)

def checkd(p):
    return os.path.isdir(p)

def mkdir(d):
    os.makedirs(d)

def getfn(p):
    return os.path.basename(p)
def getd(p):
    dirn = os.path.dirname(p)

    return dirn if dirn else "."
# def touch(p):
    # open(p, "w").close()
def get_rm_prefix(p): # a.b.c.d return a.b.c
    return ".".join(getfn(p).split(".")[0:-1])

def get_lm_prefix(p): # a.b.c.d return a
    return getfn(p).split(".")[0]



def assess_pb(man, ref, fofn, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_pb):
    if skip_pb == 1:
        return 
    jobs = []
    out_fns = []
    with open(fofn, "r") as f:
        for fl in f:
            fl_strip = fl.strip()
            fn_prefix = get_lm_prefix(getfn(fl_strip))
            out_fn = "{0}/{1}.paf".format(out_dir, fn_prefix)
            out_fns.append(out_fn)
            
            jcmd = "minimap2 -x map-pb -t {0} {1} {2} >{3}".format(str(core_lim), ref, fl_strip, out_fn)
            jjn = "minimap_{}".format(fn_prefix)
            jout = "{0}/{1}_%J.o".format(out_dir, jjn)
            jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
            
            j = hpc("lsf", cmd=jcmd, core=core_lim, mem = mem_lim, jn=jjn, out=jout, err=jerr)
            jobs.append(j)
        f.close()
    rtn = man.start(jobs)
    
    jobs = []
    if not rtn: 
        jcmd = "{0}/ast_pb -O {2} {1} >{2}/pb.bed".format(bin_dir, " ".join(out_fns), out_dir)
        jjn = "ast_pb_{}".format(spid) 
        jout = "{0}/{1}.o".format(out_dir, jjn)
        jerr = "{0}/{1}.e".format(out_dir, jjn)

        j = hpc("lsf", cmd=jcmd, core=1, mem = mem_lim, jn=jjn, out=jout, err=jerr)
        jobs.append(j)
        man.start(jobs, True)


def bwa_index(p):
    [man, ref, out_dir, spid] = p
    jobs = []
    jcmd = "bwa index {}".format(ref)
    mem_lim = 20000
    jjn = "bwa_index_{}".format(spid) 
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)

    j = hpc("lsf", cmd=jcmd, core=1, mem = mem_lim, jn=jjn, out=jout, err=jerr)
    
    jobs.append(j)
    jcmd = "bwa index {}/split.fa".format(out_dir)
    mem_lim = 20000
    jjn = "bwa_index_split_{}".format(spid) 
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)

    j = hpc("lsf", cmd=jcmd, core=1, mem = mem_lim, jn=jjn, out=jout, err=jerr)
    jobs.append(j)
    return man.start(jobs)
    # print ("code {}".format(code)) 
    # rtn.append(man.start(jobs))
    # rtn.append(code)

def aln_10x(p):
    [man, ref, fns, core_lim, mem_lim, out_fn] = p
    out_dir = getd(out_fn)
    
    [r1, r2] = fns.split('\t')
    
    pref = getfn(r1).split(".")[0]
    

    # jcmd = " ".join(['10x', '-p', pref, '-o', out_dir, '-c', r1, r2])
    # jjn = "10x_trim"
    # jout = "{0}/{1}.o".format(out_dir, jjn)
    # jerr = "{0}/{1}.e".format(out_dir, jjn)
    # j = hpc("lsf", cmd=jcmd, jn=jjn, out=jout, err=jerr)
    # rtn = man.start([j])
    
    jcmd = "bwa mem -t {0} {1} {2} {3} | samtools view -b -o - >{4}".format(core_lim, ref, r1, r2, out_fn)
     
    jjn = "bwa_mem_{}".format(pref)
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
    
    j = hpc("lsf", cmd=jcmd, core=core_lim, mem = mem_lim, jn=jjn, out=jout, err=jerr)
    rtn = man.start([j])

    return rtn
    

def trim_aln_10x(p):
    [man, ref, fns, core_lim, mem_lim, queue, out_fn, skip_trim] = p
    out_dir = getd(out_fn)
    
    [r1, r2] = fns.split('\t')
    
    pref = getfn(r1).split(".")[0]
    
    jobs = []
    rtn = 0
    if skip_trim != 1:
        jcmd = " ".join(['10x', '-p', pref, '-o', out_dir, '-c', r1, r2])
        jjn = "10x_trim"
        jout = "{0}/{1}.o".format(out_dir, jjn)
        jerr = "{0}/{1}.e".format(out_dir, jjn)
        j = hpc("lsf", cmd=jcmd, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j])
    
    if not rtn:
        if skip_trim == 1:
            jcmd = "bwa mem -t {0} {1} {2} {3} | samtools view -b -o - >{4}".format(core_lim, ref, r1, r2, out_fn)
        else:
            jcmd = "bwa mem -t {0} {1} {2}/{3}_1.fq.gz {2}/{3}_2.fq.gz | samtools view -b -o - >{4}".format(core_lim, ref, out_dir, pref, out_fn)
         
        jjn = "bwa_mem_{}".format(pref)
        jout = "{0}/{1}_%J.o".format(out_dir, jjn)
        jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
        
        j = hpc("lsf", cmd=jcmd, core=core_lim, mem = mem_lim, queue=queue, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j])

    return rtn


def assess_10x(man, ref, fofn, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_trim, skip_10x):
    # trim and align 10x
    if skip_10x == 1:
        return 
    out_fns = []
    in_fns = []
    params = []
    with open(fofn, "r") as f:
        for fl in f:
            fl_strip = fl.strip()
            # out_fn = "{0}/{1}.bam".format(out_dir, getfn(fl_strip).split(".")[0])
            out_fn = "{0}/{1}.bam".format(out_dir, get_lm_prefix(getfn(fl_strip.split('\t')[0])))
            out_fns.append(out_fn)
            in_fns.append(fl_strip)
        f.close()
  
    for i in range(len(in_fns)):
        params.append([man, ref, in_fns[i], core_lim, mem_lim, queue, out_fns[i], skip_trim])
    
    procs = []
    pl = Pool(processes=len(in_fns))
    rtvs = pl.map(trim_aln_10x, params)

    
    suc = all(v == 0 for v in rtvs)
    if suc:
        
        jcmd = "{0}/ast_10x -O {2} {2}/gaps.bed {1} >{2}/10x.bed".format(bin_dir, " ".join(out_fns), out_dir)
         
        jjn = "ast_10x_{}".format(spid)
        jout = "{0}/{1}_%J.o".format(out_dir, jjn)
        jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
        
        j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, queue="long", jn=jjn, out=jout, err=jerr)
         
        man.start([j], True)
        # man.start([j])


def assess_hic2(man, ref, fofn, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_hic):
    if skip_hic == 1:
        return 
    jobs = []
    out_fns = []

    split_ref = "{}/split.fa".format(out_dir)
    with open(fofn) as f:
        for fl in f:
            [r1, r2] = fl.strip().split('\t')
            pref = get_lm_prefix(getfn(r1))
            out_fn = "{0}/{1}.split.bam".format(out_dir, pref)
            out_fns.append(out_fn)

            jcmd = "bwa mem -SP -B10 -t {0} {1} {2} {3} | samtools view -b -o - >{4}".format(core_lim, split_ref, r1, r2, out_fn)
            jjn = "bwa_mem_{}".format(pref)
            jout = "{0}/{1}_%J.o".format(out_dir, jjn)
            jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
            j = hpc("lsf", cmd=jcmd, core=core_lim, mem = mem_lim, queue=queue, jn=jjn, out=jout, err=jerr)
            jobs.append(j)

        f.close()
    if not man.start(jobs):
        jcmd = ['samtools', 'faidx', split_ref]
        jjn = "faidx"
        jout = "{0}/{1}.o".format(out_dir, jjn)
        jerr = "{0}/{1}.e".format(out_dir, jjn)
        j = hpc(cmd=jcmd, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j])
        if rtn:
            return 
        faidx_fn = "{}.fai".format(split_ref)
         
        jcmd = "{0}/col_conts {1} >{2}/links.mat".format(bin_dir, " ".join(out_fns), out_dir)
         
        jjn = "col_conts_{}".format(spid)
        jout = "{0}/{1}_%J.o".format(out_dir, jjn)
        jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
        
        j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, jn=jjn, out=jout, err=jerr)
        
        rtn = man.start([j])
        if rtn:
            return 
        jcmd = "{0}/ast_hic2 {3} {2}/links.mat >{2}/hic2.bed".format(bin_dir, " ".join(out_fns), out_dir, faidx_fn)
         
        jjn = "ast_hic2_{}".format(spid)
        jout = "{0}/{1}_%J.o".format(out_dir, jjn)
        jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
        
        j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, jn=jjn, out=jout, err=jerr)
         
        man.start([j], True)

def assess_hic(man, ref, fofn, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_hic):
    if skip_hic == 1:
        return 
    jobs = []
    out_fns = []

    with open(fofn) as f:
        for fl in f:
            [r1, r2] = fl.strip().split('\t')
            pref = get_lm_prefix(getfn(r1))
            out_fn = "{0}/{1}.bam".format(out_dir, pref)
            out_fns.append(out_fn)

            jcmd = "bwa mem -SP -B10 -t {0} {1} {2} {3} | samtools view -b -o - >{4}".format(core_lim, ref, r1, r2, out_fn)
            jjn = "bwa_mem_{}".format(pref)
            jout = "{0}/{1}_%J.o".format(out_dir, jjn)
            jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
            j = hpc("lsf", cmd=jcmd, core=core_lim, mem = mem_lim, queue=queue, jn=jjn, out=jout, err=jerr)
            jobs.append(j)

        f.close()
    if not man.start(jobs):
        jcmd = ['samtools', 'faidx', ref]
        jjn = "faidx"
        jout = "{0}/{1}.o".format(out_dir, jjn)
        jerr = "{0}/{1}.e".format(out_dir, jjn)
        j = hpc(cmd=jcmd, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j])
        if rtn:
            return 
        faidx_fn = "{}.fai".format(ref)
        
        jcmd = "{0}/ast_hic -O {2} {2}/gaps.bed {1} >{2}/hic.bed".format(bin_dir, " ".join(out_fns), out_dir)
         
        jjn = "ast_hic_{}".format(spid)
        jout = "{0}/{1}_%J.o".format(out_dir, jjn)
        jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
        
        j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, jn=jjn, out=jout, err=jerr)
         
        man.start([j], True)

def assess_10x_hic(man, ref, fofn, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_trim, skip_10x, skip_hic):
    # p = Process(target=bwa_index, args=(man, ref, out_dir, spid, rtn)) 
    pl = Pool(processes=1)
    rtvs = pl.map(bwa_index, [[man,ref, out_dir, spid]])
    procs = []
    if rtvs[0] == 0:
        p = Process(target=assess_10x, args=(man, ref, fofn[0], core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_trim, skip_10x))
        procs.append(p)
        # p = Process(target=assess_hic, args=(man, ref, fofn[1], core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_hic))
        # procs.append(p)
        p = Process(target=assess_hic2, args=(man, ref, fofn[1], core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_hic))
        procs.append(p)
        for p in procs:
            p.start()
        for p in procs:
            p.join()


def get_bn_details(p):
    fn = getfn(p)
    fn_list = fn.split('_')
    return [fn[0], fn_list[1], fn_list[2]]


def assess_bnx_core(p):
    [man, mem_lim, core_lim, queue, ref, fn, out_dir, out_fn, bin_dir, spid, ishap] = p
    [fn_pref, tech, enzyme] = get_bn_details(get_rm_prefix(fn))
    ind = i

    # solve_dir = "/nfs/users/nfs_d/dg30/luster_dg30/dg30/projects/vgp/tools/Solve3.4_06042019a/"
    # jcmd = " ".join(["perl", solve_dir+"HybridScaffold/06042019/scripts/fa2cmap_multi_color.pl", "-e", enzyme, "1", "-i", ref, '-o', out_dir+"/fa2cmap"])
    # solve_dir = "/nfs/users/nfs_d/dg30/luster_dg30/dg30/projects/vgp/tools/Solve3.2.1_04122018/"
    # jcmd = " ".join(["perl", solve_dir+"HybridScaffold/04122018/scripts/fa2cmap_multi_color.pl", "-e", enzyme, "1", "-i", ref, '-o', out_dir+"/fa2cmap"])
    solve_dir = "/nfs/users/nfs_d/dg30/luster_dg30/dg30/projects/vgp/tools/Solve3.3_10252018/"
    jcmd = " ".join(["perl", solve_dir+"HybridScaffold/10252018/scripts/fa2cmap_multi_color.pl", "-e", enzyme, "1", "-i", ref, '-o', out_dir+"/fa2cmap"])

    jjn = "fa2cmap"
    jout = "{0}/{1}.o".format(out_dir, ind)
    jerr = "{0}/{1}.e".format(out_dir, ind)
    
    j = hpc("lsf", cmd=jcmd, mem=1000, jn=jjn, out=jout, err=jerr)
    rtn = man.start([j])
    print ("fa2cmap return value {}".format(rtn))
    if not rtn:
        jcmd = ['cp', fn, out_dir] 
        j = hpc(cmd=jcmd, out="{}/cp.o".format(out_dir))
        rtn = man.start([j])
        if not rtn:
            ref_prefix = get_rm_prefix(ref)
            ref_cmap = "{0}/fa2cmap/{1}_{2}_0kb_0labels.cmap".format(out_dir, ref_prefix, enzyme.upper())
            key_fn = "{0}/fa2cmap/{1}_{2}_0kb_0labels_key.txt".format(out_dir, ref_prefix, enzyme.upper())
            query_map = "{0}/{1}".format(out_dir, getfn(fn))
            optn = "DLE1_{}".format(tech.lower()) if enzyme == "DLE1" else tech.lower()
            # jcmd = " ".join(["python2",  solve_dir+"Pipeline/06042019/align_bnx_to_cmap.py","--prefix", enzyme, "--mol", query_map,  "--ref", ref_cmap, "--ra", solve_dir+"RefAligner/8949.9232rel/", "--nthreads", "12", "--pipeline", solve_dir+"Pipeline/06042019/", "--optArgs", solve_dir+"RefAligner/8949.9232rel/optArguments_haplotype_{}.xml".format(optn), "--output", out_dir + "/alignref_{}".format(enzyme)])
            # jcmd = " ".join(["python2",  solve_dir+"Pipeline/10252018/align_bnx_to_cmap.py","--prefix", enzyme, "--mol", query_map,  "--ref", ref_cmap, "--ra", solve_dir+"RefAligner/7437.7523rel/", "--nthreads", "12", "--pipeline", solve_dir+"Pipeline/04122018/", "--optArgs", solve_dir+"RefAligner/7437.7523rel/optArguments_haplotype_{}.xml".format(optn), "--output", out_dir + "/alignref_{}".format(enzyme)])
            jcmd = " ".join(["python2",  solve_dir+"Pipeline/10252018/align_bnx_to_cmap.py","--prefix", enzyme, "--mol", query_map,  "--ref", ref_cmap, "--ra", solve_dir+"RefAligner/7915.7989rel/", "--nthreads", "12", "--pipeline", solve_dir+"Pipeline/10252018/", "--optArgs", solve_dir+"RefAligner/7915.7989rel/optArguments_{1}_{0}.xml".format(optn, "haplotype" if ishap else "nonhaplotype"), "--output", out_dir + "/alignref_{}".format(enzyme)])
            j = hpc("lsf", cpu="avx", mem=mem_lim, core=12, queue=queue, cmd=jcmd, jn="bnx_refalign", out="{0}/bnx_refalign_{1}.o".format(out_dir, enzyme[0:4]))
            rtn = man.start([j]) 
            if not rtn:
                ref_lm_pref = get_lm_prefix(ref)
                map_path = "{0}/alignref_{2}/{1}".format(out_dir, "contigs/alignmolvref/merge/exp_refineFinal1",  enzyme)
                rmap_fn = "{}_r.cmap".format(map_path)
                qmap_fn = "{}_q.cmap".format(map_path)
                xmap_fn = "{}.xmap".format(map_path)

                jcmd = " ".join([bin_dir+"/ast_bion_bnx", rmap_fn, qmap_fn, xmap_fn, key_fn, '-O', out_dir + "/alignref_{}".format(enzyme)])
                jcmd += " >"+out_fn + " 2>" + "{0}/{1}_{2}.bed".format(out_dir, tech, enzyme)
                jn = "ast_bion_bnx"
                jout = "{0}/{1}.o".format(out_dir, jn) 
                jerr = "{0}/{1}.e".format(out_dir, jn)
                j = hpc("lsf", mem=5000, cmd=jcmd, jn="ast_bion_bnx", out=jout, err=jerr)
                rtn = man.start([j], True) 
    return rtn
def assess_bnx(man, ref, fofn_list, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_bn, ishap):
    if skip_bn == 1:
        return 
    procs = []
    out_fns = []
    params = [] 
    for i in range(len(fofn_list)):
        out_fn = "{0}/mbionano_{1}.bed".format(out_dir, i)
        out_fns.append(out_fn)
        params.append([man, mem_lim, core_lim, queue, ref, fofn_list[i], out_dir, out_fn, bin_dir, spid, ishap])
    
    pl = Pool(processes=len(fofn_list))
    rtvn = pl.map(assess_bnx_core, params)
    suc = all(v==0 for v in rtvn)
    if suc:
        if len(out_fns) > 1:
            jcmd = [bin_dir+'/union']
            jcmd.extend(out_fns)
            jout = "{}/mbn.bed".format(out_dir)
            j = hpc(cmd = jcmd, out=jout, jn="union")
            man.start([j], True)
        else:
            jcmd = ['cp', out_fns[0], out_dir+"/mbn.bed"]
            j = hpc(cmd = jcmd, jn="cp")
            man.start([j], True)

def assess_bn_core(p):
    [man, mem_lim, core_lim, queue, ref, fn, out_dir, out_fn, bin_dir, spid, ishap] = p
    [fn_pref, tech, enzyme] = get_bn_details(get_rm_prefix(fn))
    ind = i

    solve_dir = "/nfs/users/nfs_d/dg30/luster_dg30/dg30/projects/vgp/tools/Solve3.2.1_04122018/"
    jcmd = " ".join(["perl", solve_dir+"HybridScaffold/04122018/scripts/fa2cmap.pl", "-n", enzyme[0:4], "-i", ref, '-o', out_dir])
    
    jjn = "fa2cmap"
    jout = "{0}/{1}.o".format(out_dir, ind)
    jerr = "{0}/{1}.e".format(out_dir, ind)
    
    j = hpc("lsf", cmd=jcmd, mem=1000, jn=jjn, out=jout, err=jerr)
    rtn = man.start([j], True)
    if not rtn:
        jcmd = ['cp', fn, out_dir] 
        j = hpc(cmd=jcmd, out="{}/cp.o".format(out_dir))
        rtn = man.start([j])
        if not rtn:
            ref_prefix = get_rm_prefix(ref)
            ref_cmap = "{0}/fa2cmap/{1}_{2}_0Kb_0labels.cmap".format(out_dir, ref_prefix, enzyme)
            key_fn = "{0}/fa2cmap/{1}_{2}_0Kb_0labels_key.txt".format(out_dir, ref_prefix, enzyme)
            query_cmap = "{0}/{1}".format(out_dir, getfn(fn))
            optn = "DLE1_{}".format(tech.lower()) if enzyme == "DLE1" else tech.lower()
            jcmd = " ".join(["python2",  solve_dir+"Pipeline/04122018/runCharacterize.py","-t",  solve_dir+"RefAligner/7437.7523rel/RefAligner","-q", query_cmap, "-r", ref_cmap, "-p", solve_dir+"Pipeline/04122018/", "-a", solve_dir+"RefAligner/7437.7523rel/optArguments_{1}_{0}.xml".format(optn, "haplotype" if ishap else "nonhaplotype"), "-n","2"])
             
            j = hpc("lsf", cpu="avx", mem=mem_lim, core=core_lim, queue=queue, cmd=jcmd, jn="bn_refalign", out="{0}/bn_refalign_{1}.o".format(out_dir, enzyme[0:4]))
            rtn = man.start([j], True) 
            if not rtn:
                ref_lm_pref = get_lm_prefix(ref)
                map_path = "{0}/alignref/{1}".format(out_dir, get_rm_prefix(fn))
                rmap_fn = "{}_r.cmap".format(map_path)
                qmap_fn = "{}_q.cmap".format(map_path)
                xmap_fn = "{}.xmap".format(map_path)

                jcmd = " ".join([bin_dir+"/ast_bion", rmap_fn, qmap_fn, xmap_fn, key_fn])
                jcmd += " >"+out_fn + " 2>" + "{0}/{1}_{2}.bed".format(out_dir, tech, enzyme)
                jn = "ast_bion"
                jout = "{0}/{1}.o".format(out_dir, jn) 
                jerr = "{0}/{1}.e".format(out_dir, jn)
                j = hpc("lsf", mem=5000, cmd=jcmd, jn="ast_bion", out=jout, err=jerr)
                rtn = man.start([j], True) 
    return rtn
def assess_bn(man, ref, fofn_list, core_lim, mem_lim, queue, out_dir, bin_dir, spid, skip_bn, ishap):
    if skip_bn == 1:
        return 
    procs = []
    out_fns = []
    params = [] 
    for i in range(len(fofn_list)):
        out_fn = "{0}/bionano_{1}.bed".format(out_dir, i)
        out_fns.append(out_fn)
        params.append([man, mem_lim, core_lim, queue, ref, fofn_list[i], out_dir, out_fn, bin_dir, spid, ishap])
    
    pl = Pool(processes=len(fofn_list))
    rtvn = pl.map(assess_bn_core, params)
    suc = all(v==0 for v in rtvn)
    if suc:
        if len(out_fns) > 1:
            jcmd = [bin_dir+'/union']
            jcmd.extend(out_fns)
            jout = "{}/bn.bed".format(out_dir)
            j = hpc(cmd = jcmd, out=jout, jn="union")
            man.start([j], True)
        else:
            jcmd = ['cp', out_fns[0], out_dir+"/bn.bed"]
            j = hpc(cmd = jcmd, jn="cp")
            man.start([j], True)
            
def acc(man, ref, out_dir, bin_dir, spid):
    # if not checkf():
        # jcmd = ['samtools', 'faidx', ref]
        # jjn = "faidx"
        # jout = "{0}/{1}.o".format(out_dir, jjn)
        # jerr = "{0}/{1}.e".format(out_dir, jjn)
        # j = hpc(cmd=jcmd, jn=jjn, out=jout, err=jerr)
        # rtn = man.start([j])
        # if rtn:
            # return 
    beds = []
    for fn in ["gaps.bed", "10x.bed", "bn.bed", "hic.bed", "pb.bed"]:
        fpath = "{0}/{1}".format(out_dir,fn)
        if checkf(fpath):
            beds.append(fpath)
    if len(beds):
        jcmd = "{0}/acc {1} > {2}/acc.bed".format(bin_dir, " ".join(beds), out_dir) 
        jjn = "acc_{}".format(spid)
        jout = "{0}/acc_{1}.o".format(out_dir, spid)
        jerr = "{0}/acc_{1}.e".format(out_dir, spid)
        j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, jn=jjn, out=jout, err=jerr)
        man.start([j], True)
    #acc contig
    
    jcmd = "{0}/acc {2}/gaps.bed {2}/pb.bed {2}/bn.bed > {2}/pb_bn.bed".format(bin_dir, " ".join(beds), out_dir) 
    jjn = "acc_contig_{}".format(spid)
    jout = "{0}/acc_contig_{1}.o".format(out_dir, spid)
    jerr = "{0}/acc_contig_{1}.e".format(out_dir, spid)
    j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, jn=jjn, out=jout, err=jerr)
    man.start([j], True)

    jcmd = "{0}/acc {2}/gaps.bed {2}/10x.bed {2}/hic2.bed {2}/bn.bed > {2}/10x_hic2_bn.bed".format(bin_dir, " ".join(beds), out_dir) 
    jjn = "acc_scaf_{}".format(spid)
    jout = "{0}/acc_scaf_{1}.o".format(out_dir, spid)
    jerr = "{0}/acc_scaf_{1}.e".format(out_dir, spid)
    j = hpc("lsf", cmd=jcmd, core=1, mem = 20000, jn=jjn, out=jout, err=jerr)
    man.start([j], True)
    


def punchlist(man, ref, out_dir, bin_dir, spid):
    # all_bed = "{}/acc.bed".format(out_dir)
    # if checkf(all_bed):
        # jcmd = "{0}/pchlst {1}/gaps.bed {2} > {1}/pchlst.bed".format(bin_dir, out_dir, all_bed) 
        # jjn = "pchlst_{}".format(spid)
        # jout = "{0}/pchlst_{1}.o".format(out_dir, spid)
        # jerr = "{0}/pchlst_{1}.e".format(out_dir, spid)
        # j = hpc(cmd = jcmd, out=jout, err=jerr, jn=jjn)
        # man.start([j], True)
    contig_acc = "{}/pb_bn.bed".format(out_dir)
    if checkf(contig_acc):
        jcmd = "{0}/pchlst -c {1}/gaps.bed {2} > {1}/pchlst_ctg.bed".format(bin_dir, out_dir, contig_acc) 
        jjn = "pchlst_ctg_{}".format(spid)
        jout = "{0}/pchlst_ctg_{1}.o".format(out_dir, spid)
        jerr = "{0}/pchlst_ctg_{1}.e".format(out_dir, spid)
        j = hpc(cmd = jcmd, out=jout, err=jerr, jn=jjn)
        man.start([j], True)
    scaf_acc = "{}/10x_hic2_bn.bed".format(out_dir)
    if checkf(scaf_acc):
        jcmd = "{0}/pchlst  {1}/gaps.bed {2} > {1}/pchlst_scaf.bed".format(bin_dir, out_dir, scaf_acc) 
        jjn = "pchlst_scf_{}".format(spid)
        jout = "{0}/pchlst_scf_{1}.o".format(out_dir, spid)
        jerr = "{0}/pchlst_scf_{1}.e".format(out_dir, spid)
        j = hpc(cmd = jcmd, out=jout, err=jerr, jn=jjn)
        man.start([j], True)

    ctg_pchlst = "{}/pchlst_ctg.bed".format(out_dir)
    scf_pchlst = "{}/pchlst_scaf.bed".format(out_dir)
    if checkf(scaf_acc) and checkf(contig_acc):
        faidx_fn = "{}.fai".format(ref)
        jcmd = "{0}/union_brks -x {1}/gaps.bed {2} {3} > {1}/pchlst_chrom.bed".format(bin_dir, out_dir, ctg_pchlst, scf_pchlst) 
        jjn = "pchlst_chrom_{}".format(spid)
        jout = "{0}/pchlst_chrom_{1}.o".format(out_dir, spid)
        jerr = "{0}/pchlst_chrom_{1}.e".format(out_dir, spid)
        j = hpc(cmd = jcmd, out=jout, err=jerr, jn=jjn)
        man.start([j], True)

def postproc(man, fofn, out_dir, bin_dir, spid):
    pchlst_bed = "{}/pchlst.bed".format(out_dir)
    in_fns = []
    with open(fofn, "r") as f:
        for fl in f:
            fl_strip = fl.strip()
            fn_prefix = get_lm_prefix(getfn(fl_strip)) 
            in_fn = "{0}/{1}.paf".format(out_dir, fn_prefix)
            in_fns.append(in_fn)
        f.close()

    if checkf(pchlst_bed):
        jcmd = "{0}/ast_postproc {1} {2} > {3}/post_pchlst.bed".format(bin_dir, pchlst_bed, " ".join(in_fns), out_dir)
        jout = "{}/postproc.o".format(out_dir)
        jerr = "{}/postproc.e".format(out_dir)
        j = hpc("lsf", cmd=jcmd, core=1, mem = 3000, jn="ast_postproc", out=jout, err=jerr)
        man.start([j])

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print ("run <config> <bins> <id>")
        sys.exit(1)
    else:
        config_fn = sys.argv[1]
        bin_dir = sys.argv[2]
        spid = sys.argv[3]

        f = open(config_fn, "r")
        config_dict = json.load(f)
        
        out_dir = config_dict["out_dir"]
        ref = config_dict["ref"]
        if not checkd(out_dir):
            mkdir(out_dir)


        man = manager(retries=2) 
        
        jcmd = "{0}/detgaps {1} > {2}/gaps.bed".format(bin_dir, ref, out_dir)
        jout = "{}/detgaps.o".format(out_dir)
        jerr = "{}/detgaps.e".format(out_dir)
        
        j = hpc("lsf", cmd=jcmd, jn="detgaps", out=jout, err=jerr)
        if man.start([j], True):
            print ("fail to generate gaps for {}".format(ref))
            sys.exit(1)
        jcmd = "{0}/split_fa {1} > {2}/split.fa".format(bin_dir, ref, out_dir)
        jout = "{}/split_ref.o".format(out_dir)
        jerr = "{}/split_ref.e".format(out_dir)
        
        j = hpc("lsf", cmd=jcmd, jn="split_ref", out=jout, err=jerr)
        if man.start([j], True):
            print ("fail to split_ref for {}".format(ref))
            sys.exit(1)

        procs = []
        
        # func_list = [assess_pb, assess_10x_hic, assess_bn, assess_bnx]
        func_list = [assess_pb, assess_10x_hic, assess_bn]
        # func_list = [aassess_bn]
        # key_list = ["pb", "10x_hic", "bn", "bnx"]
        key_list = ["pb", "10x_hic", "bn"]
        # key_list = ["bn"]
        for i in range(len(func_list)):
            if key_list[i] not in config_dict: 
                continue
            cur_d = config_dict[key_list[i]]
            if i == 1:
                p = Process(target=func_list[i], args=(man, ref, cur_d["fofn"], cur_d["core"], cur_d["mem"], cur_d["queue"], out_dir, bin_dir, spid, cur_d["skip_trim"] if "skip_trim" in cur_d else 0, cur_d["skip_10x"] if "skip_10x" in cur_d else 0, cur_d["skip_hic"] if "skip_hic" in cur_d else 0))
            elif i == 0:
                skip_str = "skip_" + key_list[i]
                p = Process(target=func_list[i], args=(man, ref, cur_d["fofn"], cur_d["core"], cur_d["mem"], cur_d["queue"], out_dir, bin_dir, spid, cur_d[skip_str] if skip_str in cur_d else 0))
            else:
                skip_str = "skip_" + key_list[i]
                p = Process(target=func_list[i], args=(man, ref, cur_d["fofn"], cur_d["core"], cur_d["mem"], cur_d["queue"], out_dir, bin_dir, spid, cur_d[skip_str] if skip_str in cur_d else 0, cur_d["nonhap"] if "nonhap" in cur_d else 1))
            procs.append(p)
        for p in procs:
            p.start()
        for p in procs:
            p.join()
        # gen  accumulate file
        p = Process(target=acc, args=(man, ref, out_dir, bin_dir, spid))

        p.start()
        
        p.join()
        # gen punch list
        p = Process(target=punchlist, args=(man, ref, out_dir, bin_dir, spid))

        p.start()
        
        p.join()
        
        # p = Process(target=postproc, args=(man, config_dict["pb"]["fofn"], out_dir, bin_dir, spid))
        # p.start()
        
        # p.join()
