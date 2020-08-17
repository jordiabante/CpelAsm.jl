# CpelAsm
![CI](https://github.com/jordiabante/CpelAsm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelAsm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelAsm.jl/blob/master/LICENSE.md)

## Generate GFF files

The first step consists in generating the 2 necessary GFF files: heterozygous and homozygous.
The former is the one that contains the haplotypes analyzed, while the latter contains the
homozygous regions of the genome.

In the following example, the maximum number of CpG sites allowed per haplotype is 25 (`n_max=25`).
In addition, CpelAsm extends the window delimited by the first and last SNP in the haplotype by 100
bp left and right (`win_exp=100`).

```julia
# Parameters
n_max = 25
win_exp = 100

# Paths
project_path = "/path/to/SampleName"
vcf = "$(project_path)/vcf/SampleName.phased.vcf"
fasta = "$(project_path)/fasta/SampleName.masked_hg19.fa"
het_gff = "$(project_path)/CpelAsmOut/SampleName_het.cpelasm.gff"
hom_gff = "$(project_path)/CpelAsmOut/SampleName_hom.cpelasm.gff"

# Call
gen_gffs([het_gff,hom_gff],fasta,vcf,win_exp,n_max)

```

## BedGraphs MML, NME, & UC

The next step consists in estimating the allele-specific Ising models in the haplotypes
and generating bedGraph files with MML1/2, NME1/2, and UC. CpelAsm can be parallelized 
when multiple CPUs are available by loading first the package `Distributed` and then 
loading CpelAsm through the macro ``@everywhere`.

In the following example, the maximum size of a subregion is fixed to 500 (`g_max=500`), 
the minimum average depth is set to 8 (`cov_ths=8`), and the WGBS reads are trimmed 5 bp
on each end (`trim=(5,5,5,5)`).

```julia
# Deps
using Distributed
@everywhere using CpelAsm

# Parameters
g_max = 500
cov_ths = 8
trim = (5,5,5,5)

# Paths
project_path = "/path/to/SampleName"
fasta = "$(project_path)/fasta/SampleName.masked_hg19.fa"
het_gff = "$(project_path)/cpelasm/SampleName_het.cpelasm.gff"
bam1 = "$(project_path)/bam/SampleName.sort.genome1.bam"
bam2 = "$(project_path)/bam/SampleName.sort.genome2.bam"
uc_path = "$(project_path)/cpelasm/SampleName_uc.bedGraph"
mml1_path = "$(project_path)/cpelasm/SampleName_mml1.bedGraph"
mml2_path = "$(project_path)/cpelasm/SampleName_mml2.bedGraph"
nme1_path = "$(project_path)/cpelasm/SampleName_nme1.bedGraph"
nme2_path = "$(project_path)/cpelasm/SampleName_nme2.bedGraph"
tobs_path = [mml1_path,mml2_path,nme1_path,nme2_path,uc_path]

# Call
comp_tobs(bam1,bam2,het_gff,fasta,tobs_path;g_max=g_max,cov_ths=cov_ths,trim=trim)

```

## Generate Null Statistics

The next step consists in generating null statistics to be able to perform hypothesis
testing. As discussed in the previous point, CpelAsm can be parallelized when multiple 
CPUs are available by first loading the package `Distributed` and then loading CpelAsm 
through the macro ``@everywhere`.

In the following example, the maximum size of a subregion is fixed to 500 (`g_max=500`), 
the minimum average depth is set to 8 (`cov_ths=8`), and the WGBS reads are trimmed 5 bp
on each end (`trim=(5,5,5,5)`).

```julia
# Deps
using Distributed
@everywhere using CpelAsm

# Parameters
cov_ths = 8
g_max = 500
trim = (5,5,5,5)

# Paths
project_path = "/path/to/SampleName"
fasta = "$(project_path)/fasta/SampleName.masked_hg19.fa"
bam1 = "$(project_path)/bam/SampleName.sort.genome1.bam"
bam2 = "$(project_path)/bam/SampleName.sort.genome2.bam"
het_gff = "$(project_path)/cpelasm/SampleName_het.cpelasm.gff"
hom_gff = "$(project_path)/cpelasm/SampleName_hom.cpelasm.gff"
null_tpdm_path = "$(project_path)/cpelasm/SampleName_tpdm_null.bedGraph"
null_tmml_path = "$(project_path)/cpelasm/SampleName_tmml_null.bedGraph"
null_tnme_path = "$(project_path)/cpelasm/SampleName_tnme_null.bedGraph"
tnull_path = [null_tmml_path,null_tnme_path,null_tpdm_path]
uc_path = "$(project_path)/cpelasm/SampleName_uc.bedGraph"
mml1_path = "$(project_path)/cpelasm/SampleName_mml1.bedGraph"
mml2_path = "$(project_path)/cpelasm/SampleName_mml2.bedGraph"
nme1_path = "$(project_path)/cpelasm/SampleName_nme1.bedGraph"
nme2_path = "$(project_path)/cpelasm/SampleName_nme2.bedGraph"
tobs_path = [mml1_path,mml2_path,nme1_path,nme2_path,uc_path]

# Call
comp_tnull(bam,het_gff,hom_gff,fasta,tobs_path,tnull_path;
    g_max=g_max,cov_ths=cov_ths,trim=trim,n_max=n_max)

```

## Perform Allele-Specific Methylation Detection

The final step is to compute a p-value for each haplotype for each one of the three
quantities (Tmml, Tnme, Tpdm). The returned p-values are corrected using the
Benjamini-Hochberg and allow for control of the false discovery rate (FDR).

In the following example, hypothesis testing is performed for haplotypes with at
most 25 CpG sites (`n_max=25`), and a minimum of 1,000 null statistics to perform
hypothesis testing.

```julia
# Deps
using CpelAsm

# Parameters
n_max = 25
n_null = 1000

# Paths
project_path = "/path/to/SampleName"
fasta = "$(project_path)/fasta/SampleName.masked_hg19.fa"
bam1 = "$(project_path)/bam/SampleName.sort.genome1.bam"
bam2 = "$(project_path)/bam/SampleName.sort.genome2.bam"
null_tpdm_path = "$(project_path)/cpelasm/SampleName_tpdm_null.bedGraph"
null_tmml_path = "$(project_path)/cpelasm/SampleName_tmml_null.bedGraph"
null_tnme_path = "$(project_path)/cpelasm/SampleName_tnme_null.bedGraph"
tnull_path = [null_tmml_path,null_tnme_path,null_tuc_path]
uc_path = "$(project_path)/cpelasm/SampleName_uc.bedGraph"
mml1_path = "$(project_path)/cpelasm/SampleName_mml1.bedGraph"
mml2_path = "$(project_path)/cpelasm/SampleName_mml2.bedGraph"
nme1_path = "$(project_path)/cpelasm/SampleName_nme1.bedGraph"
nme2_path = "$(project_path)/cpelasm/SampleName_nme2.bedGraph"
tobs_path = [mml1_path,mml2_path,nme1_path,nme2_path,uc_path]
pVal_tmml_path = "$(project_path)/cpelasm/SampleName_tmml_pvals.bedGraph"
pVal_tnme_path = "$(project_path)/cpelasm/SampleName_tnme_pvals.bedGraph"
pVal_tpdm_path = "$(project_path)/cpelasm/SampleName_tpdm_pvals.bedGraph"
pVal_path = [pVal_tmml_path,pVal_tnme_path,pVal_tuc_path]

# Call
comp_pvals(tobs_path,tnull_path,p_path,n_max,n_null)

```

## Running CpelAsm with a single command

The following example shows how to perform all the steps shown above in a single command.

```julia
# Deps
using Distributed
@everywhere using CpelAsm

# Parameters
cov_ths = 8
g_max = 500
win_exp = 100
trim = (5,5,5,5)

# Paths
project_path = "/path/to/SampleName"
vcf = "$(project_path)/vcf/SampleName.phased.vcf"
fasta = "$(project_path)/fasta/SampleName.masked_hg19.fa"
bam1 = "$(project_path)/bam/SampleName.sort.genome1.bam"
bam2 = "$(project_path)/bam/SampleName.sort.genome2.bam"
bamu = "$(project_path)/bam/SampleName.sort.unassigned.bam"
outdir = "$(project_path)/cpelasm/"

# Call
run_analysis(bam1,bam2,bamu,vcf,fasta,outdir;g_max=g_max,cov_ths=cov_ths,
    win_exp=win_exp,trim=trim,n_null=n_null,n_max=n_max)

```
