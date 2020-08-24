# CpelAsm
![CI](https://github.com/jordiabante/CpelAsm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelAsm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelAsm.jl/blob/master/LICENSE.md)

## Toy Example

The package includes a small toy example, which should take no more
than 20 seconds to run, for illustrative purposes. The example in
particular consists of 2 haplotypes in a single chromosome of an
artificial reference genome. Each haplotype has two alleles: `a1` and
`a2`. For simplicity, CpG sites in both alleles are distributed
independently according to a Bernoulli distribution. In `a1`, the
probability of methylation is `p1=0.8`, while in `a2` the probability of
methylation is `p2=0.2`. Thus, CpG sites in the former have a mean
methylation level (MML) of 0.8, while that of the later is 0.2. Given the
symmetry of the problem, however, both alleles have the same Shannon
entropy (see Shannon entropy of a flip coin). Thus differential analysis
only identifies ASM imbalances in for statistics Tmml and Tpdm. That is,
CpelAsm only identifies mean methylation level differences and probability
distribution of methylation differences. The output bedGraph files with
the suffix `pvals` contain the results of the statistical test performed
for each haplotype. These files can be found in the `out/` directory. To
run the toy example run the following commands in a `julia`'s REPL:

```julia
# Load CpelAsm
using CpelAsm

# Define I/O
dir = "/path/to/CpelAsm.jl/test/"
b1 = "$(dir)/bam/example.a1.bam"
b2 = "$(dir)/bam/example.a2.bam"
fa = "$(dir)/fasta/n-masked/example.fa"
vcf = "$(dir)/vcf/example.vcf"
out = "$(dir)/out/"

# Run CpelAsm
run_analysis(b1,b2,b1,vcf,fa,out;g_max=50,win_exp=10,n_null=50,n_max=10)
```
