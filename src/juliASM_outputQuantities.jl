"""
    `comp_mml(N,α,β)`

Function that computes mean methylation level assuming the Ising model for
a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> comp_mml(4,0.0,0.0)
0.5
```
"""
function comp_mml(N::Int64,a::Float64, b::Float64)::Float64
    # Generate Xcal
    xcal = generate_xcal(N)

    # Traverse xcal
    mml = 0.0
    for x in xcal
        mml += sum(0.5*(x.+1)) * comp_lkhd(x, a, b)
    end

    # Return
    return 1.0/N * mml
end # end comp_mml
"""
    `comp_shanH(N,α,β)`

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
    return -1.0/log2(N+1) * shanH
end # end comp_shanH
"""
    `comp_mi(POS_CG,η1,η2)`

Function that computes the mutual information between haplotype and methylation
state assuming an Ising model for each allele with parameters η1 and η2
respectively. POS_CG contains the homozygous CpG sites, as well as the CpG sites
for each haplotype.

# Examples
```julia-repl
julia> comp_mi([[1],[1],[1]],[-10.0,-10.0],[10.0,10.0])
0.9999999375540618
```
"""
function comp_mi(cpg_pos::Array{Array{Int64,1},1},t1::Array{Float64,1},
                 t2::Array{Float64,1})::Float64

    # Generate Xcal
    xcal = generate_xcal(length(cpg_pos[1]))
    hom_in_h1 = findall(x->x in cpg_pos[1], cpg_pos[2])
    hom_in_h2 = findall(x->x in cpg_pos[1], cpg_pos[3])

    # Traverse xcal
    mi = 0.0
    lkhd1 = 0.0
    lkhd2 = 0.0
    for x in xcal
        x1 = zeros(Int64,length(cpg_pos[2]))
        x2 = zeros(Int64,length(cpg_pos[3]))
        x1[hom_in_h1] = x
        x2[hom_in_h2] = x
        lkhd1 = comp_lkhd(x1, t1[1], t1[2])
        lkhd2 = comp_lkhd(x2, t2[1], t2[2])
        ave_lkhd = 0.5 * (lkhd1+lkhd2)
        mi += lkhd1 * log2(lkhd1/ave_lkhd) + lkhd2 * log2(lkhd2/ave_lkhd)
    end

    # Return
    return 0.5 * mi
end # end comp_mi
