# JuliASM

A `julia` package for allele-specific methylation (ASM) analysis.

## Description

`JuliASM` is based on the method in [CITE], and it draws ideas from
statistical physics and information theory to identify allele-specific
methylation events at the haplotype level.

## Getting Started

### Prerequisites

* `julia v1.0.0` or greater.
* `git`.

### Installing

`JuliASM` and dependencies can be installed via the following command:
```julia
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/jordiabante/JuliASM.jl.git"))'
```

## Running the tests

In a `julia` session run
```julia
(v1.0) pkg> test JuliASM
```

## Toy Example

The package includes a small toy example for illustrative purposes.
The example consists of two alleles `a1` and `a2`. The former has a
mean-methylation level (MML) of 0.4, while that of the later is 0.6.
Nevertheless, given the symmetry of the problem, both alleles have
the same Shannon entropy. Thus differential  analysis only identifies
differences in terms of MML. The output bedGraph files can be found
in `out_path`. To run the toy example run the following commands in
a `julia` session:

```julia
using JuliASM
dir = "/path/to/JuliASM.jl/test/"
b1 = "$(dir)/bam/example.a1.bam"
b2 = "$(dir)/bam/example.a2.bam"
fa = "$(dir)/fasta/example.fa"
vcf = "$(dir)/vcf/example.vcf"
out = "$(dir)/out/"
run_analysis(b1,b2,b1,vcf,fa,out;blk_size=20,cov_ths=3,cov_b=20.0,win_exp=10,mc_null=1000,n_max=5)
```

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md)
file for details.

## References
[CITE]
