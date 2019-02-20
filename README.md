# JuliASM

A `julia` package for allele-specific methylation analysis.

## Description

`JuliASM` is based on the method in XXX, and it draws ideas from statistical
physics and information theory to identify significant differences between
the methylation patterns of each haplotype.

## Getting Started

### Prerequisites

* `julia v1.0.0` or greater.
* `git`.

### Installing

`JuliASM` and dependencies can be installed via the following command:
```julia
julia -e 'using Pkg; Pkg.clone("git@github.com:jordiabante/JuliASM.jl.git")'
```

## Running the tests

In a `julia` session run
```julia
pkg> test JuliASM
```

## Toy Example

The package includes a small toy example for illustrative purposes. The
example consists of two alleles `a1` and `a2`. The former has a
mean-methylation level of 0.9, while the later has 0.1. The mutual
information is almost equal to 1 in all variants, and the permutation test
returns p-values≈0. The output bedGraph files can be found in `out_path`.
To run the toy example run the following commands in a `julia` session:

```julia
using JuliASM
dir = "/path/to/JuliASM.jl/test/"
bam1_path="$(dir)/bam/example.a1.bam"
bam2_path="$(dir)/bam/example.a2.bam"
fasta_path = "$(dir)/fasta/example.fa"
vcf_path = "$(dir)/vcf/example.vcf"
out_path="$(dir)/out/"
run_asm_analysis(bam1_path,bam2_path,vcf_path,fasta_path,out_path)
```

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md)
file for details.

## Acknowledgments

* Inspiration
