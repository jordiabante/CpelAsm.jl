"""
    comp_mml(N,α,β)

Function that computes mean methylation level assuming the Ising model for
a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> comp_mml(4,0.0,0.0)
0.5
```
"""
function comp_mml(N::Int64,a::Float64, b::Float64)::Float64
    # TODO: find more efficient way to compute it.
    # Generate Xcal
    xcal = generate_xcal(N)

    # Traverse xcal
    lkhd = 0.0
    mml = 0.0
    for x in xcal
        lkhd = comp_lkhd(x, a, b)
        mml += sum(0.5*(x.+1)) * lkhd
    end

    # Return
    return 1/N * mml
end # end comp_mml
"""
    comp_shanH(N,α,β)

Function that computes Shannon's entropy (in bits) assuming the Ising model for
a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> comp_shanH(4,0.0,0.0)
4.0
```
"""
function comp_shanH(N::Int64,a::Float64, b::Float64)::Float64
    # Generate Xcal
    xcal = generate_xcal(N)

    # Traverse xcal
    lkhd = 0.0
    shanH = 0.0
    for x in xcal
        lkhd = comp_lkhd(x, a, b)
        shanH += lkhd * log2(lkhd)
    end

    # Return
    return -1/log2(N+1) * shanH
end # end comp_shanH
"""
    comp_mi(N,θ1,θ2)

Function that computes the mutual informtion between the allele of origin and
the methylation state vector of size N assuming an Ising model for each allele
with parameters θ1 and θ2 respectively.

# Examples
```julia-repl
julia> comp_mi(4,[-10.0,-10.0],[10.0,10.0])
1.0000000000000382
```
"""
function comp_mi(N::Int64,t1::Array{Float64,1},t2::Array{Float64,1})::Float64
    # Generate Xcal
    xcal = generate_xcal(N)

    # Traverse xcal
    lkhd1 = 0.0
    lkhd2 = 0.0
    mi = 0.0
    for x in xcal
        lkhd1 = comp_lkhd(x, t1[1], t1[2])
        lkhd2 = comp_lkhd(x, t2[1], t2[2])
        ave_lkhd = 0.5 * (lkhd1+lkhd2)
        mi += lkhd1 * log2(lkhd1/ave_lkhd) + lkhd2 * log2(lkhd2/ave_lkhd)
    end

    # Return
    return 0.5 * mi
end # end comp_mi
