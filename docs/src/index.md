# CpelAsm

[![Build Status](https://travis-ci.com/jordiabante/CpelAsm.jl.svg?token=XZfbD5CqoU7r1YJmmbNE&branch=master)](https://travis-ci.com/jordiabante/CpelAsm.jl)

## Description

CpelAsm is a julia package for haplotype allele-specific methylation, based
on the method in [CITE], that draws ideas from statistical physics and
information theory to identify allele-specific methylation events at the
haplotype level.

## Testing

CpelAsm is tested against Julia `1.X` on Linux and OS X.

## Installation

### Prerequisites

* julia > v1.0.0
* git.

### Command

`CpelAsm` and dependencies can be installed using the following command:
```julia
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/jordiabante/CpelAsm.jl.git"))'
```

### Installation tests

In a `julia` session run
```julia
(v1.0) pkg> test CpelAsm
```

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](../../LICENSE.md)
file for details.

## References
[CITE]
