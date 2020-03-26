# CpelAsm

![CI](https://github.com/jordiabante/CpelAsm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelAsm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelAsm.jl/blob/master/LICENSE.md)

## Description

CpelAsm is a julia package especifically desgined for haplotype allele-specific 
methylation based on the method in [1]. CpelAsm draws ideas from statistical 
physics and information theory to detect allele-specific methylation imbalances 
at the haplotype level.

## Testing

CpelAsm is tested against Julia `1.3.0` on the latest versions of Linux, macOS and Windows.

## Installation

### Prerequisites

* julia v1.3.0
* git.

### Command

`CpelAsm` and dependencies can be installed using the following command:
```julia
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/jordiabante/CpelAsm.jl.git"))'
```

### Installation tests

In a `julia` session run
```julia
(v1.3) pkg> test CpelAsm
```

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](../../LICENSE.md)
file for details.

## References
[1] Abante, J., Fang, Y., Feinberg, A.P., Goutsias, J., Detection of haplotype-dependent 
allele-speciﬁc DNA methylation in WGBS data, *Nature Communications* 2020 XYZ.