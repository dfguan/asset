<hr>
<h1 align="center", style="color:red;">Bullitin</h1>
<hr>
<ol>
<li>UI modifications</li>
<table>
  <tr>
    <th>Program</th>
    <th>Operation</th>
    <th>Old</th>
    <th>New</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>ast_bion_bnx</td>
    <td>Update</td>
    <td>m</td>
    <td>c</td>
    <td>minimum coverage</td>
  </tr>
    <tr>
    <td>ast_bion_bnx</td>
    <td>Update</td>
    <td>M</td>
    <td>C</td>
    <td>maximum coverage</td>
  </tr>
   <td>ast_pb</td>
    <td>Update</td>
    <td>m</td>
    <td>c</td>
    <td>minimum coverage</td>
  </tr>
   <td>ast_pb</td>
    <td>Update</td>
    <td>M</td>
    <td>C</td>
    <td>minimum coverage</td>
  </tr>
    </tr>
   <td>ast_pb</td>
    <td>Add</td>
    <td></td>
    <td>r</td>
    <td>minimum coverage ratio</td>
  </tr>
      </tr>
   <td>ast_pb</td>
    <td>Add</td>
    <td></td>
    <td>R</td>
    <td>maximum coverage ratio</td>
  </tr>
    </tr>
   <td>ast_pb</td>
    <td>Add</td>
    <td></td>
    <td>s</td>
    <td>use soft limit</td>
  </tr>
</table>
<li>New functions</li>
<ul>
<li>ast_pb supports soft limits now. the lower coverage bound for a contig/scaffold is calculated by max(min_ratio * ave, min_cov), and the upper bound is min(max_ratio * ave, max_cov), where <b>ave</b> is the average coverage for the contig/scaffold. 
</ul>
</ol>
<hr>

# Asset
Assembly evaluation tool

# Workflow
Asset can use sequencing data from four platforms (Pacbio, 10X, Bionano, HiC) to accumulate the support evidence for a *de novo* assembly. You can use [Pipeline Guide](#pg) to build your own pipeline. 


# Dependencies
## C Library
- zlib

## Third Party Tools
- minimap2
- bwa 
- samtools
- RefAligner
- runner (optional, only necessary when you run my pipeline script run_asset.py)


# Installation

```
git clone https://github.com/dfguan/asset.git
cd asset/src && make
```
There will be bin direcotory keeping all asset executables under asset after compiling successfully.

# <a name="pg"> Pipeline Guide </a>
## Preprocessing
Given an assembly *asm*, use the following command to preprocess.

```
bin/detgaps $asm > $output_dir/gaps.bed
``` 

## Pacbio Processing
Given a pacbio files list *pblist* and the assembly *asm*, apply the following command to get Pacbio support regions.

```
for fl in $pblist
do
	minimap2 -xmap-pb -t 12 $asm $fl > $fl.paf
done

bin/ast_pb $fl1.paf $fl2.paf $fl3.paf ... >$output_dir/pb.bed 2>ast_pb.log
```

## 10X Processing
Given a 10x read files list *10xlist* (suppose in fastq.gz format) and the assembly *asm*, use the following command to get 10X support regions.

```
bwa index $asm
while read -r r1 r2
do
	prefix=`basename $r1 .fastq.gz`
	bin/10x -c -p $prefix $r1 $r2 # generate trimmed read files $prefix_{1,2}.fq.gz
	bwa mem -t12 $asm $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -b - > $prefix.bam
done < $10xlist

bin/ast_10x $output_dir/gaps.bed $bam1 $bam2 $bam3 ... >$output_dir/10x.bed 2>ast_10x.log

```
**10x_trim** is available at [dfguan/utls](https://github.com/dfguan/utls).
## Bionano Processing
### Consensus map (.cmap)
Given a bionano consensus map files list *bnlist* (suppose in .cmap format) and the assembly *asm*, use the following command to get Bionano support regions.

```

solve_dir=/nfs/users/nfs_d/dg30/luster_dg30/dg30/projects/vgp/tools/Solve3.2.1_04122018/
for fl in $bnlist
do
	fn=`basename $fl`
	fn_pref=`echo $fn | cut -d_ -f1`
	tech=`echo $fn | cut -d_ -f2`
	enzyme=`echo $fn | cut -d_ -f3`
	perl $solve_dir/HybridScaffold/04122018/scripts/fa2cmap.pl -n ${enzyme:0:4} -i ref -o $output_dir
	cp $fl $output_dir
	ref_prefix=${asm%.*}
	ref_cmap=$output_dir/fa2cmap/"$ref_prefix"_"$enzyme"_0Kb_0labels.cmap
	key_fn=$output_dir/fa2cmap/"$ref_prefix"_"$enzyme"_0Kb_0labels_key.txt
	query_cmap=$output_dir/$fn
	optn=${tech,,}
	if [ "$enzyme" = "DLE1" ]
	then
		optn="DLE1_"$optn
	fi
   python2 $solve_dir/Pipeline/04122018/runCharacterize.py -t   $solve_dir/RefAligner/7437.7523rel/RefAligner -q $query_cmap -r  $ref_cmap -p $solve_dir/Pipeline/04122018/ -a $solve_dir/RefAligner/7437.7523rel/optArguments_nonhaplotype_"$optn".xml -n 4
	map_path=$output_dir/alignref/${fn%.*} 
	rmap_fn="$map_path"_r.cmap
	qmap_fn="$map_path"_q.cmap
	xmap_fn="$map_path".xmap
	bin/ast_bion $rmap_fn $qmap_fn $xmap_fn $key_fn > $output_dir/bionano_"$tech"_"$enzyme".bed 2>ast_bion_"$tech"_"$enzyme".log
done
fln=`wc -l $bnlist | awk {print $1}`
# only consider one or two bionano files
if [ "$fln" -eq 2 ]
then
	bin/union bionano_*.bed > bn.bed
else
	cp bionano_*.bed bn.bed
fi
```

### Molecular map (.bnx)

Update soon

## HiC Processing 
Given a HiC files list *hiclist* (suppose in fastq.gz format) and the assembly *asm*, use the following command to get Bionano support regions.

```
bin/split_fa $asm > split.fa
samtools faidx split.fa 
bwa index split.fa
while read -r r1 r2
do
	prefix=`basename $r1 .fastq.gz`
	bwa mem -SP -B10 -t12 split.fa $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -b - > $prefix.bam
done < $hiclist
bin/col_conts *.bam > $output_dir/links.mat
bin/ast_hic2 split.fa.fai $output_dir/links.mat >$output_dir/hic2.bed 2>ast_hic.log

```

## Evidence Accumulation 
once you got the bed files from Pacbio (pb.bed), 10X (10x.bed), Bionano (bn.bed), HiC (hic.bed), run the following command to get accumulation bed file and detect break points.

```
bin/acc $output_dir/gaps.bed $output_dir/{pb,bn}.bed $output_dir/bn.bed > $output_dir/pb_bn.bed 
bin/acc $output_dir/gaps.bed $output_dir/{10x,hic2,bn}.bed > $output_dir/10x_hic2_bn.bed  
```

## Mis-assemblies Detection
```
bin/pchlst -c $output_dir/gaps.bed $output_dir/pb_bn.bed > $output_dir/pchlst_ctg.bed
bin/pchlst $output_dir/gaps.bed $output_dir/10x_hic2_bn.bed > $output_dir/pchlst_scaf.bed 
bin/union_brks $output_dir/gaps.bed $output_dir/pchlst_{ctg,scaf}.bed > $output_dir/pchlst_final.bed
```

## Contact
jWellcome to use, you can use github webpage to report an issue or email me dfguan9@gmail.com with any advice.
