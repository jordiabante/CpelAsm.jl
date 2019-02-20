###############################################################################
# CONSTANTS
###############################################################################
const LOG2 = log(2)                      # Ln(2)
###############################################################################
# FUNCTIONS
###############################################################################
"""
`comp_ex([N1,...,NK],[α1,...,αK],β)`

Function that computes the mean methylation vector E[X] assuming the Ising model
for a methylation state vector of size [N1,...,NK] and parameters [α1,...,αK]
and β.

# Examples
```julia-repl
julia> JuliASM.comp_ex([4],[0.0],0.0)
4-element Array{Float64,1}:
0.5
0.5
0.5
0.5
```
"""
function comp_ex(n::Array{Int64,1},a::Array{Float64,1},b::Float64)::Array{Float64,1}

    # Loop over all positions
    x = zeros(Int64,sum(n))
    y = zeros(Float64,sum(n))
    @inbounds for j in 1:sum(n)
        x[j] = 1
        y[j] = comp_lkhd(x,n,a,b)
        x[j] = 0
    end

    # Return
    return 2*y.-1

end # end comp_ex
"""
`comp_exx([N1,...,NK],[α1,...,αK],β;DIST)`

Function that computes mean methylation vector E[X_{i}X_{+LAG}}] assuming the
Ising model for a methylation state vector of size [N1,...,NK] and parameters
[α1,...,αK] and β.

# Examples
```julia-repl
julia> JuliASM.comp_exx([4],[0.0],0.0)
3-element Array{Float64,1}:
 0.0
 2.220446049250313e-16
 0.0
```
"""
function comp_exx(n::Array{Int64,1},a::Array{Float64,1},b::Float64;dist::Int64=1)::Array{Float64,1}

    # Loop over all positions
    x = zeros(Int64,sum(n))
    y = zeros(Float64,sum(n)-dist)
    @inbounds for j in 1:(sum(n)-dist)
        x[j] = x[j+dist] = 1
        y[j] = comp_lkhd(x,n,a,b) + comp_lkhd(-x,n,a,b)
        x[j] = x[j+dist] = 0
    end

    # Return
    return 2*y.-1

end # end comp_exx
"""
    `comp_cov([N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that returns the covariance matrix of a methylation vector given the
[N1,...,NK], [α1,...,αK], β, E[X], and E[XX].

# Examples
```julia-repl
julia> JuliASM.comp_cov([4],[0.0],0.0,ex,exx)
??
```
"""
function comp_cov(n::Array{Int64,1},a::Array{Float64,1},b::Float64,ex::Array{Float64,1},exx::Array{Float64,1})::Array{Float64,2}

    # Initialize matrix
    ntot = sum(n)
    cov = zeros(Float64,ntot,ntot)

    # Add first subdiagonal E[X_{i}X_{i+1}]
    cov[2:(ntot+1):((ntot-1)*ntot)] = exx

    # Add covariance terms E[X_{i}X_{i+k-1}]
    @inbounds for k=3:ntot
        cov[k:(ntot+1):((ntot-k+1)*ntot)] = comp_exx(n,a,b;dist=k-1)
    end

    # Symmetrize
    cov += transpose(cov)

    # Substract E[X_i]E[X_j]
    cov -= ex * ex'

    # Add 1's to diagonal of matrix
    cov[1:(ntot+1):ntot^2] += ones(Float64,ntot)

    # Return Σ
    return cov

end # end comp_cov
"""
    `comp_evec(Σ)`

Function that returns the eigenvector associated with the largest eigenvalue of
the covariance matrix Σ.

# Examples
```julia-repl
julia> JuliASM.comp_evec(cov)
??
```
"""
function comp_evec(cov::Array{Float64,2})::Array{Float64,1}

    # Return evec
    return abs.(eigvecs(cov)[:,size(cov)[1]])

end # end comp_evec
"""
    `comp_mml(EX)`

Function that computes mean methylation level per CpG site assuming the Ising
model for a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> JuliASM.comp_mml(JuliASM.comp_ex([4],[0.0],0.0))
0.5
```
"""
function comp_mml(ex::Array{Float64,1})::Float64

    # Return
    return abs(round(0.5/length(ex)*sum(ex)+0.5;digits=4))

end # end comp_mml
"""
    `comp_corr(EX,EXX)`

Function that returns the correlation vector between consecutive CpG sites given
the vector of E[X] and E[XX].

# Examples
```julia-repl
julia> JuliASM.comp_corr(JuliASM.comp_ex([4],[0.0],0.0),JuliASM.comp_exx([4],[0.0],0.0))
3-element Array{Float64,1}:
 0.0
 8.881784197001252e-16
 0.0
```
"""
function comp_corr(ex::Array{Float64,1},exx::Array{Float64,1})::Array{Float64,1}

    # Initialize vector of 1's
    ov = ones(Float64,length(ex)-1)

    # Return ρ
    return (exx.-ex[1:(end-1)].*ex[2:end]) ./ sqrt.((ov.-ex[1:(end-1)].^2).*(ov.-ex[2:end].^2))

end # end comp_corr
"""
    `comp_nme([N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that computes normalized methylation entropy assuming the Ising model for
a methylation state vector of size [N1,...,NK] and parameters [α1,...,αK] and β.

# Examples
```julia-repl
julia> N=[10]
julia> a=[0.0]
julia> b=0.0
julia> JuliASM.comp_nme(N,a,b,JuliASM.comp_ex(N,a,b),JuliASM.comp_exx(N,a,b))
1.0
```
"""
function comp_nme(n::Array{Int64,1},a::Array{Float64,1},b::Float64,ex::Array{Float64,1},exx::Array{Float64,1})::Float64

    # Vector of α
    avec = vcat([a[i]*ones(Float64,n[i]) for i in 1:length(n)]...)

    # Return
    return abs(round((log(comp_Z(n,a,b))-avec'*ex-b*sum(exx))/(sum(n)*LOG2);digits=4))

end # end comp_nme
"""
    `comp_nmi(nme0,nme1,nme2)`

Function that computes the normalized mutual information between haplotype and
methylation state assuming an Ising model for both single-allele and allele-
agnostic models.

# Examples
```julia-repl
julia> JuliASM.comp_nmi(1.0,0.5,1.5)
0.0
```
"""
function comp_nmi(h::Float64,h1::Float64,h2::Float64)::Float64

    # Return
    return abs(round(h - 0.5*(h1+h2);digits=4))

end # end comp_nmi
"""
    `bifurcate(xpool,sel)`

Function that divides xpool into two mutually exclusive sets based on sel.

# Examples
```julia-repl
julia> JuliASM.bifurcate([[1,1],[1,-1],[-1,1]], [1,2])
(Array{Int64,1}[[1, 1], [1, -1]], Array{Int64,1}[[-1, 1]])
```
"""
function bifurcate(xpool::Array{Array{Int64,1},1},sel::Vector{T}) where T <: Integer
    x = xpool[sel]
    asel = trues(length(xpool))
    asel[sel] .= false
    y = xpool[asel]
    return x,y
end
"""
    `perm_test(XOBS1,XOBS2,TOBS,POS_CG)`

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs1=JuliASM.gen_ising_full_data(100,4;a=1.0);
julia> xobs2=JuliASM.gen_ising_full_data(100,4;a=-1.0);
julia> eta1=JuliASM.est_eta(xobs1);
julia> eta2=JuliASM.est_eta(xobs2);
julia> hom=[1,5,10,15];
julia> cpg_pos=[hom,hom,hom];
julia> tobs=JuliASM.comp_mi(cpg_pos,eta1,eta2);
julia> JuliASM.perm_test(xobs1,xobs2,tobs,cpg_pos)
0.0
```
"""
function perm_test(xobs1::Array{Array{Int64,1},1},xobs2::Array{Array{Int64,1},1},
                   tobs::Float64,cpg_pos::Array{Array{Int64,1},1})::Float64
    # Initialize
    better = worse = 0.0
    xpool = vcat(xobs1, xobs2)

    # Loop over combinations and compute statistic
    i = 1
    for subset in combinations(1:length(xpool), length(xobs2))
      # Bifurcate
      test, control = bifurcate(xpool, subset)

      # Estimate parameters
      length(cpg_pos[2])==1 ? eta1=est_alpha(control) : eta1=est_eta(control)
      length(cpg_pos[3])==1 ? eta2=est_alpha(test) : eta2=est_eta(test)

      # Compute MI for partition and compared to that observed
      comp_mi(cpg_pos,eta1,eta2)>tobs ? better+=1.0 : worse+=1.0

      # If enough permutations leave
      i<100 ? i+=1 : break
    end

    # Return p-value
    return better/(worse+better)
end # end perm_test
