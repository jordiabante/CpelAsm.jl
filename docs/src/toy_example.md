# CpelAsm

[![Build Status](https://travis-ci.com/jordiabante/CpelAsm.jl.svg?token=XZfbD5CqoU7r1YJmmbNE&branch=master)](https://travis-ci.com/jordiabante/CpelAsm.jl)

## Toy Example

The package includes a small toy example for illustrative purposes.
The example consists of two alleles `a1` and `a2`. The former has a
mean-methylation level (MML) of 0.4, while that of the later is 0.6.
Nevertheless, given the symmetry of the problem, both alleles have
the same Shannon entropy. Thus differential  analysis only identifies
differences in terms of MML. The output bedGraph files can be found
in `out_path`. To run the toy example run the following commands in
a `julia`'s REPL:

```julia
using CpelAsm
dir = "/path/to/CpelAsm.jl/test/"
b1 = "$(dir)/bam/example.a1.bam"
b2 = "$(dir)/bam/example.a2.bam"
fa = "$(dir)/fasta/example.fa"
vcf = "$(dir)/vcf/example.vcf"
out = "$(dir)/out/"
run_analysis(b1,b2,b1,vcf,fa,out;g_max=25,cov_ths=3,cov_b=20.0,win_exp=10,n_null=100,n_max=10)
```
