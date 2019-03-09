###################################################################################################
# CONSTANTS
###################################################################################################
const LOG2 = log(2)                         # Ln(2)
const XCAL1 = [-1,1]                        # 𝒳 with N=1
const XCAL2 = [[-1,-1],[-1,1],[1,-1],[1,1]] # 𝒳 with N=2
###################################################################################################
# USED FUNCTIONS
###################################################################################################
"""
`comp_ex([N1,...,NK],[α1,...,αK],β)`

Function that computes the mean methylation vector E[X] assuming the Ising model for a methylation
state vector of size `[N1,...,NK]` and parameters `[α1,...,αK]` and `β`.

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
`comp_exx([N1,...,NK],[α1,...,αK],β;r=1)`

Function that computes mean methylation vector E[X_{i}X_{i+r}}] assuming the Ising model for a
methylation state vector of size `[N1,...,NK]` and parameters `[α1,...,αK]` and `β`.

# Examples
```julia-repl
julia> JuliASM.comp_exx([4],[0.0],0.0)
3-element Array{Float64,1}:
 0.0
 2.220446049250313e-16
 0.0
```
"""
function comp_exx(n::Array{Int64,1},a::Array{Float64,1},b::Float64;r::Int64=1)::Array{Float64,1}

    # Loop over all positions
    x = zeros(Int64,sum(n))
    y = zeros(Float64,sum(n)-r)
    @inbounds for j in 1:(sum(n)-r)
        x[j] = x[j+r] = 1
        y[j] = comp_lkhd(x,n,a,b) + comp_lkhd(-x,n,a,b)
        x[j] = x[j+r] = 0
    end

    # Return
    return 2*y.-1

end # end comp_exx
"""
    `comp_mml(EX)`

Function that computes mean methylation level (MML) given the first order moment E[X].

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
    `comp_exlng(Z,[N1,...,NK],[α1,...,αK],β)`

Function that returns the third summand in the normalized methylation entropy (NME) defined
for the homozygous CpG sites, with positions determined by binary vector Z, assuming an
allele-specific vector with [N1,...,NK] CpG sites, with parameters [α1,...,αK] and β.

# Examples
```julia-repl
julia> n=[10]
julia> z=trues(10);z[5]=false
julia> a=[0.0]
julia> b=0.0
julia> JuliASM.comp_exlng(z,n,a,b)
0.6931471805599452
```
"""
function comp_exlng(z::BitArray{1},n::Vector{Int64},a::Vector{Float64},b::Float64)::Float64

    # Find changes to/from 0 in x vector
    hetst = findall(isequal(-1),abs.(z[2:end]) - abs.(z[1:(end-1)])) .+ 1
    hetend = findall(isequal(1),abs.(z[2:end]) - abs.(z[1:(end-1)])) .+ 1

    # Determine whether it starts/finishes with true/false
    z[1]==false && pushfirst!(hetst,1)
    z[end]==false && push!(hetend,sum(n)+1)

    # Find subregion label for each CpG site
    subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)

    # Loop over heterozygous stretches
    y = 0.0
    @inbounds for i=1:length(hetst)

        # Choose 𝒳 depending on boundary or not
        bound = hetst[i]==1 || hetend[i]==sum(n)+1 ? true : false
        xcal = bound ? XCAL1 : XCAL2

        # Figure out b (block IDs of heterozygous) and r (α indices)
        bid = subid[hetst[i]:(hetend[i]-1)]
        nhet = [count(x->x==id,bid) for id in unique(bid)]

        # Compute expectation
        for x in xcal

            # Find αp1 and αp2
            xaug = zeros(Int64,sum(n))
            if bound
                # If heterozygous CpG on boundary of X
                hetst[i]==1 ? xaug[hetend[i]]+=x : xaug[hetst[i]-1]+=x
                ap1 = hetst[i]==1 ? a[1] : 2.0*x*b+a[subid[hetst[i]]]
                ap2 = hetend[i]==sum(n)+1 ? a[end] : 2.0*x*b+a[subid[hetend[i]-1]]
            else
                # If heterozygous CpG not on boundary of X
                xaug[[hetst[i]-1,hetend[i]]]+=x
                ap1 = 2.0*x[1]*b+a[subid[hetst[i]]]
                ap2 = 2.0*x[2]*b+a[subid[hetend[i]-1]]
            end

            # Add contribution to total sum of expectations
            y += log(comp_g(nhet,a[unique(bid)],b,ap1,ap2))*comp_lkhd(xaug,n,a,b)

        end

    end

    # Return sum of all expectations
    return y

end # end comp_exlng
"""
    `comp_nme(Z,[N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that computes normalized methylation entropy (NME) over the CpG sites determined by binary
vector Z. This is done by assuming an Ising model for the allele-specific methylation state vector
of size `[N1,...,NK]`, parameters `[α1,...,αK]` and `β`.

# Examples
```julia-repl
julia> n=[10]
julia> z=trues(10);z[5]=false;
julia> a=[0.0]
julia> b=0.0
julia> JuliASM.comp_nme(z,n,a,b,JuliASM.comp_ex(n,a,b),JuliASM.comp_exx(n,a,b))
1.0
```
"""
function comp_nme(z::BitArray{1},n::Vector{Int64},a::Vector{Float64},b::Float64,
                  ex::Vector{Float64},exx::Vector{Float64})::Float64

    # Define output quantity
    h = 0.0

    # Vector of α for homozygous
    avec = vcat([a[i]*ones(Float64,n[i]) for i in 1:length(n)]...)[z]

    # 1st term (log partition function)
    h += log(comp_Z(n,a,b))

    # 2nd term (average potential energy over homozygous CpG sites)
    h -= avec'*ex[z]+b*sum(exx[z[1:(end-1)].*z[2:end]])

    # 3rd term if necessary
    h -= all(z) ? 0.0 : comp_exlng(z,n,a,b)

    # Return nme
    return abs(round(h/(sum(z)*LOG2);digits=4))

end # end comp_nme
"""
    `comp_uc(nme0,nme1,nme2)`

Function that computes the uncertainty coefficient (UC) between haplotype and methylation
state assuming an Ising model for both single-allele and allele-agnostic models.

# Examples
```julia-repl
julia> JuliASM.comp_uc(1.0,0.5,1.5)
0.0
```
"""
function comp_uc(nme0::Float64,nme1::Float64,nme2::Float64)::Float64

    # Return
    return abs(round(1.0-0.5*(nme1+nme2)/nme0;digits=4))

end # end comp_uc
###################################################################################################
# UNUSED FUNCTIONS
###################################################################################################
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

      # TODO: Compute MI for partition and compared to that observed
      comp_mi(cpg_pos,eta1,eta2)>tobs ? better+=1.0 : worse+=1.0

      # If enough permutations leave
      i<100 ? i+=1 : break
    end

    # Return p-value
    return better/(worse+better)

end # end perm_test
"""
    `comp_cov([N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that returns the covariance matrix of a methylation vector given the `[N1,...,NK]`,
`[α1,...,αK]`, `β`, E[X], and E[XX].

# Examples
```julia-repl
julia> ex = JuliASM.comp_ex([4],[0.0],0.0)
julia> exx = JuliASM.comp_exx([4],[0.0],0.0)
julia> JuliASM.comp_cov([4],[0.0],0.0,ex,exx)
4×4 Array{Float64,2}:
 1.0          0.0          2.22045e-16  0.0
 0.0          1.0          2.22045e-16  2.22045e-16
 2.22045e-16  2.22045e-16  1.0          0.0
 0.0          2.22045e-16  0.0          1.0
```
"""
function comp_cov(n::Array{Int64,1},a::Array{Float64,1},b::Float64,ex::Array{Float64,1},
                  exx::Array{Float64,1})::Array{Float64,2}

    # Initialize matrix
    ntot = sum(n)
    cov = zeros(Float64,ntot,ntot)

    # Add first subdiagonal E[X_{i}X_{i+1}]
    cov[2:(ntot+1):((ntot-1)*ntot)] = exx

    # Add covariance terms E[X_{i}X_{i+k-1}]
    @inbounds for k=3:ntot
        cov[k:(ntot+1):((ntot-k+1)*ntot)] = comp_exx(n,a,b;r=k-1)
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

Function that returns the eigenvector associated with the largest eigenvalue of the covariance
matrix Σ.

# Examples
```julia-repl
julia> ex = JuliASM.comp_ex([4],[0.0],0.0);
julia> exx = JuliASM.comp_exx([4],[0.0],0.0);
julia> cov = JuliASM.comp_cov([4],[0.0],0.0,ex,exx);
julia> JuliASM.comp_evec(cov)
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 1.0
```
"""
function comp_evec(cov::Array{Float64,2})::Array{Float64,1}

    # Return evec
    return abs.(eigvecs(cov)[:,size(cov)[1]])

end # end comp_evec
"""
    `comp_corr(EX,EXX)`

Function that returns the correlation vector between consecutive CpG sites given the vector of
first order moments E[X] and second order moments E[XX].

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
THIS FUNCTION IS JUST FOR TESTING PURPOSES!

    `comp_nme_xcal(Z,[N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that computes normalized methylation entropy (NME) over the homozygous CpG sites. This is
done by assuming an Ising model for the allele-specific methylation state vector of size
`[N1,...,NK]`, parameters `[α1,...,αK]` and `β`. The homozygous part of the allele-specific vector
is determined by binary vector Z (i.e., via Hadamard product Z*X, where * is the Hadamard product).

# Examples
```julia-repl
julia> n=[10]
julia> z=trues(sum(n))
julia> a=[0.0]
julia> b=0.0
julia> JuliASM.comp_nme_xcal(z,n,a,b,JuliASM.comp_ex(n,a,b),JuliASM.comp_exx(n,a,b))
1.0
```
"""
function comp_nme_xcal(z::BitArray{1},n::Vector{Int64},a::Vector{Float64},b::Float64,
                       ex::Vector{Float64},exx::Vector{Float64})::Float64

    # Loop over 𝒳h
    h = 0.0
    w = zeros(Int64,sum(n))
    hom_ind = findall(y->y==true,z)
    xcal = generate_xcal(length(hom_ind))
    for x in xcal
        w[hom_ind] =+ x
        lkhd = comp_lkhd(w,n,a,b)
        h += lkhd * log(lkhd)
        w[hom_ind] =- x
    end

    # Return
    return -h/(sum(z)*LOG2)

end # end comp_nme_xcal
"""
    `comp_nmi(nme0,nme1,nme2)`

Function that computes the normalized mutual information (NMI) between haplotype and methylation
state assuming an Ising model for both single-allele and allele-agnostic models.

# Examples
```julia-repl
julia> JuliASM.comp_nmi(1.0,0.5,1.5)
0.0
```
"""
function comp_nmi(nme0::Float64,nme1::Float64,nme2::Float64)::Float64

    # Return
    return abs(round(nme0-0.5*(nme1+nme2);digits=4))

end # end comp_nmi
