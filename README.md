# CpelAsm

[![Build Status](https://travis-ci.com/jordiabante/CpelAsm.jl.svg?token=XZfbD5CqoU7r1YJmmbNE&branch=master)](https://travis-ci.com/jordiabante/CpelAsm.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://jordiabante.github.io/CpelAsm.jl/stable)

## Description

`CpelAsm` is based on the method in [CITE], and it draws ideas from
statistical physics and information theory to identify allele-specific
methylation events at the haplotype level.

## Getting Started

### Prerequisites

* `julia v1.0.0` or greater.
* `git`.

### Installing

`CpelAsm` and dependencies can be installed via the following command:
```julia
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/jordiabante/CpelAsm.jl.git"))'
```

## Running the tests

In a `julia` session run
```julia
(v1.0) pkg> test CpelAsm
```

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

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md)
file for details.

## References
[CITE]
