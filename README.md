# <img src="./docs/src/assets/logo.png" width="30%" align="right" /> CpelAsm

![CI](https://github.com/jordiabante/CpelAsm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelAsm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelAsm.jl/blob/master/LICENSE.md)

## Description

CpelAsm is a julia package especifically designed for haplotype allele-specific 
methylation based on the method in [1]. CpelAsm draws ideas from statistical 
physics and information theory to detect allele-specific methylation imbalances 
at the haplotype level.

## Testing

CpelAsm is tested against Julia `1.3.0` on the latest versions of Linux, macOS and Windows.

## Getting Started

### Prerequisites

* julia v1.3.0
* git.

### Installing

`CpelAsm` and dependencies can be installed via the following command in julia's REPL:
```julia
(v1.3) pkg> add https://github.com/jordiabante/CpelAsm.jl.git
```

## Running the tests

In a `julia` session run
```julia
(v1.3) pkg> test CpelAsm
```

## Toy Example

The package includes a small toy example for illustrative purposes.
The example consists of two alleles `a1` and `a2`. The former has a
mean-methylation level (MML) of 0.8, while that of the later is 0.2.
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
fa = "$(dir)/fasta/n-masked/example.fa"
vcf = "$(dir)/vcf/example.vcf"
out = "$(dir)/out/"
run_analysis(b1,b2,b1,vcf,fa,out;g_max=50,cov_ths=5,cov_b=2.0,win_exp=10,n_null=50,n_max=10)
```

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md)
file for details.

## References
[1] Abante, J., Fang, Y., Feinberg, A.P., Goutsias, J., Detection of haplotype-dependent 
allele-speciﬁc DNA methylation in WGBS data, *Nature Communications* 2020 XYZ.
