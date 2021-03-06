###################################################################################################
# CONSTANTS
###################################################################################################
const LOG2 = log(2)                         # Ln(2)
const XCAL1 = [-1,1]                        # 𝒳 with N=1
const XCAL2 = [[-1,-1],[-1,1],[1,-1],[1,1]] # 𝒳 with N=2
###################################################################################################
# STRUCTS
###################################################################################################
mutable struct IntBitVec <: AbstractVector{Bool}
    data::UInt64
end
Base.size(::IntBitVec) = (64,)
@inline function Base.getindex(v::IntBitVec,i::Int)
   @boundscheck checkbounds(v, i)
   v.data & (1 << (i-1)) != 0
end
###################################################################################################
# USED FUNCTIONS
###################################################################################################
"""
`comp_ex([N1,...,NK],[α1,...,αK],β)`

Function that computes the mean methylation vector E[X] assuming the Ising model for a methylation
state vector of size `[N1,...,NK]` and parameters `[α1,...,αK]` and `β`.

# Examples
```julia-repl
julia> CpelAsm.comp_ex([4],[0.0],0.0)
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
```
"""
function comp_ex(n::Vector{Int64},a::Vector{Float64},b::Float64)::Vector{Float64}

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
julia> CpelAsm.comp_exx([4],[0.0],0.0)
3-element Array{Float64,1}:
 0.0
 2.220446049250313e-16
 0.0
```
"""
function comp_exx(n::Vector{Int64},a::Vector{Float64},b::Float64;r::Int64=1)::Vector{Float64}

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
    `comp_mml(Z,EX)`

Function that computes mean methylation level (MML) given at the CpG sites marked by binary vector
Z and given first order moment E[X] vector. This function is only used when not all CpG sites are
homozygous since it can compute MML over a subset of CpG sites.

# Examples
```julia-repl
julia> ex=CpelAsm.comp_ex([4],[0.0],0.0); 
julia> CpelAsm.comp_mml(trues(4),ex)
0.5
```
"""
function comp_mml(z::BitArray{1},ex::Vector{Float64})::Float64

    # Return
    return abs(round(0.5/length(ex[z])*sum(ex[z])+0.5;digits=8))

end # end comp_mml
"""
    `comp_mml_∇([N1,...,NK],∇logZ)`

Function that computes mean methylation level (MML) over the entire haplotype. This function takes
advantage of the fact that E[S(X)]=Adot(θ), and it's much faster than `comp_mml`. This function is
only used when all CpG sites are homozygous.

# Examples
```julia-repl
julia> n=[1,1,1]; θ=[1.0,-1.0,1.0,0.0];
julia> ∇logZ = CpelAsm.get_grad_logZ(n,θ);
julia> CpelAsm.comp_mml_∇(n,∇logZ)
0.62693236
```
"""
function comp_mml_∇(n::Vector{Int64},∇logZ::Vector{Float64})::Float64

    # Return
    return abs(round(0.5*(1.0+1.0/sum(n)*sum(∇logZ[1:length(n)]));digits=8))

end # end comp_mml_∇
"""
    `comp_exlng(Z,[N1,...,NK],[α1,...,αK],β)`

Function that returns the third summand in the normalized methylation entropy (NME) defined
for the homozygous CpG sites, with positions determined by binary vector Z, assuming an
allele-specific vector with [N1,...,NK] CpG sites, with parameters [α1,...,αK] and β.

# Examples
```julia-repl
julia> n=[10]; z=trues(10); z[5]=false; α=[0.0]; β=0.0;
julia> CpelAsm.comp_exlng(z,n,α,β)
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
of size `[N1,...,NK]`, parameters `[α1,...,αK]` and `β`. This function is used when not all CpG
sites are homozygous since it can compute NME over a subset of CpG sites.

# Examples
```julia-repl
julia> n=[10]; z=trues(sum(n)); z[5]=false; α=[0.0]; β=0.0;
julia> CpelAsm.comp_nme(z,n,α,β,CpelAsm.comp_ex(n,α,β),CpelAsm.comp_exx(n,α,β))
1.0
```
"""
function comp_nme(z::BitArray{1},n::Vector{Int64},a::Vector{Float64},b::Float64,ex::Vector{Float64},
                  exx::Vector{Float64})::Float64

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
    return abs(round(h/(sum(z)*LOG2);digits=8))

end # end comp_nme
"""
    `comp_nme_∇([N1,...,NK],θ,∇logZ)`

Function that computes normalized methylation entropy (NME) using the ∇logZ information. This
function can only be used in `comp_tnull` since it does not allow computation on a subset of CpG
sites. This function is only used when all CpG sites are homozygous.

# Examples
```julia-repl
julia> n=[4]; α=[0.0]; β=0.0; θ=vcat([α,β]...); ∇logZ=CpelAsm.get_grad_logZ(n,θ);
julia> CpelAsm.comp_nme_∇(n,θ,∇logZ)
1.0
```
"""
function comp_nme_∇(n::Vector{Int64},θ::Vector{Float64},∇logZ::Vector{Float64})::Float64

    # Return
    return abs(round(1.0/(sum(n)*LOG2)*(log(comp_Z(n,θ[1:(end-1)],θ[end]))-θ'*∇logZ);digits=8))

end # end comp_nme_∇
"""
    `gen_x_mc(N,α,β)`

Function that generates a methylation vector from an Ising model with `[N1,...,NK]`, and parameters
`[α1,...,αK]` and `β`.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> CpelAsm.gen_x_mc([5],[0.0],0.0)
5-element Array{Int64,1}:
  1
  1
  1
 -1
  1
```
"""
function gen_x_mc(n::Vector{Int64},a::Vector{Float64},b::Float64)::Vector{Int64}

    # Find subregion label for each CpG site
    subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)

    # Sample first CpG site
    x = Vector{Int64}()
    p = comp_lkhd(pushfirst!(zeros(Int64,sum(n)-1),1),n,a,b)
    push!(x,sample(XCAL1,Weights([1-p,p])))
    sum(n)>1 || return x

    # Sequentially sample from Markov chain
    @inbounds for i=2:sum(n)

        # Add boundary function g() if necessary
        expaux = exp(a[subid[i]]+b*x[i-1])
        if i < sum(n)
            ap1p = 2.0*b + a[subid[i]]
            ap1m = -2.0*b + a[subid[i]]
            n_miss = [count(x->x==sr,subid[(i+1):sum(n)]) for sr in unique(subid[(i+1):sum(n)])]
            gp = comp_g(n_miss,a[unique(subid[(i+1):sum(n)])],b,ap1p,a[end])
            gm = comp_g(n_miss,a[unique(subid[(i+1):sum(n)])],b,ap1m,a[end])
            p = gp*expaux/(gp*expaux+gm/expaux)
        else
            p = expaux/(expaux+1.0/expaux)
        end

        # Add i-th CpG site
        push!(x,sample(XCAL1,Weights([1-p,p])))

    end

    # Return methylation vector
    return x

end # end gen_x_mc
"""
    `comp_nme_mix_mc(Z1,Z2,N1,N2,theta1,theta2)`

Function that estimates the entropy for the Ising mixture model using Monte Carlo. This is done
assuming an Ising model for the allele-specific methylation state vector of size `[N1,...,NK]`,
parameters `[α1,...,αK]` and `β`, for each allele. The homozygous part of each model is determined
by binary vector Z (i.e., via Hadamard product Z*X, where * is the Hadamard product).

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> α1=[0.0]; β1=0.0; α2=[1.0]; β2=1.0;
julia> n1=[10]; n2=[10]; z1=trues(sum(n1)); z2=trues(sum(n2));
julia> CpelAsm.comp_nme_mix_mc(z1,z2,n1,n2,vcat(α1,β1),vcat(α2,β2))
0.6278665195419888
```
"""
function comp_nme_mix_mc(z1::BitArray{1},z2::BitArray{1},n1::Vector{Int64},n2::Vector{Int64},
                         t1::Vector{Float64},t2::Vector{Float64};L::Int64=1000)::Float64

    # For loop to sample from mixture
    h = 0.0
    @inbounds for i=1:L

        # Sample allele
        y = sample(1:2,Weights([0.5,0.5]))

        # Sample x given allele
        x = y==1 ? gen_x_mc(n1,t1[1:(end-1)],t1[end]) : gen_x_mc(n2,t2[1:(end-1)],t2[end])
        x0 = y==1 ? x[z1] : x[z2]

        # Add contribution
        x1 = zeros(Int64,sum(n1))
        x2 = zeros(Int64,sum(n2))
        x1[z1] = x0
        x2[z2] = x0
        h += log(comp_lkhd(x1,n1,t1[1:(end-1)],t1[end])+comp_lkhd(x2,n2,t2[1:(end-1)],t2[end]))

    end

    # Return
    return min(1.0,max(0.0,1.0/sum(z1)*(1.0-h/(L*LOG2))))

end # end comp_nme_mix_mc
"""
    `comp_nme_mix_exact(Z1,Z2,N1,N2,theta1,theta2)`

Function that exactly computes the NME for the Ising mixture model. This is done by assuming an
Ising model for the allele-specific methylation state vector of size `[N1,...,NK]`, parameters
`[α1,...,αK]` and `β`. The homozygous part of the allele-specific vector is determined by binary
vector Z (i.e., via Hadamard product Z*X, where * is the Hadamard product). This is done by
traversing 𝒳h.

# Examples
```julia-repl
julia> n1=[10]; n2=[10];
julia> α1=[0.0]; β1=0.0; α2=[1.0]; β2=1.0;
julia> z1=trues(sum(n1)); z2=trues(sum(n2));
julia> CpelAsm.comp_nme_mix_exact(z1,z2,n1,n2,vcat(α1,β1),vcat(α2,β2))
0.6239489595165986
```
"""
function comp_nme_mix_exact(z1::BitArray{1},z2::BitArray{1},n1::Vector{Int64},n2::Vector{Int64},
                            t1::Vector{Float64},t2::Vector{Float64})::Float64

    # Loop over 𝒳h
    h = 0.0
    n = sum(z1)
    w1 = zeros(Int64,sum(n1))
    w2 = zeros(Int64,sum(n2))
    @inbounds for i=1:(2^n)
        z = IntBitVec(i)[1:n]
        x = -ones(Int64,n)
        x[z] .= 1
        w1[z1] =+ x
        w2[z2] =+ x
        lkhd1 = comp_lkhd(w1,n1,t1[1:(end-1)],t1[end])
        lkhd2 = comp_lkhd(w2,n2,t2[1:(end-1)],t2[end])
        h += (lkhd1+lkhd2) * log(lkhd1+lkhd2)
        w1[z1] =- x
        w2[z2] =- x
    end

    # Return
    return min(1.0,max(0.0,1.0/n*(1.0-0.5*h/LOG2)))

end # end comp_nme_mix_exact
"""
    `comp_uc(Z1,Z2,N1,N2,theta1,theta2,h1,h2)`

Function that exactly computes UC over the homozygous CpG sites. This is done by assuming an Ising
model for the allele-specific methylation state vector of size `[N1,...,NK]`, parameters
`[α1,...,αK]` and `β`. The homozygous part of the allele-specific vector is determined by binary
vector Z (i.e., via Hadamard product Z*X, where * is the Hadamard product).

# Examples
```julia-repl
julia> n1=[20]; n2=[20];
julia> z1=trues(sum(n1)); z2=trues(sum(n2));
julia> α1=[2.5]; α2=[2.5]; β1=0.0; β2=0.0;
julia> h1=CpelAsm.comp_nme(z1,n1,α1,β1,CpelAsm.comp_ex(n1,α1,β1),CpelAsm.comp_exx(n1,α1,β1));
julia> h2=CpelAsm.comp_nme(z2,n2,α2,β2,CpelAsm.comp_ex(n2,α2,β2),CpelAsm.comp_exx(n2,α2,β2));
julia> CpelAsm.comp_uc(z1,z2,n1,n2,vcat(α1,β1),vcat(α2,β2),h1,h2)
0.63304398
```
"""
function comp_uc(z1::BitArray{1},z2::BitArray{1},n1::Vector{Int64},n2::Vector{Int64},
                  t1::Vector{Float64},t2::Vector{Float64},h1::Float64,h2::Float64)::Float64

    # Compute h(x)
    h = sum(z1)<17 ? comp_nme_mix_exact(z1,z2,n1,n2,t1,t2) : comp_nme_mix_mc(z1,z2,n1,n2,t1,t2)
        # print_log("h: $(h)")

    # If H(X)≈0 then UC=0
    uc = isapprox(h,0.0,atol=1e-5) ? 0.0 : min(1.0,max(0.0,1.0-0.5*(h1+h2)/h))

    # Return
    return round(uc;digits=8)

end # end comp_uc
###################################################################################################
# DEV FUNCTIONS
###################################################################################################
"""
    `comp_mmlv([N1,...,NK],∇logZ)`

Function that computes the mean methylation level vector (MMLV) for each subregion. Reporting
output for subregion might be more informative than just an overall mean. This function takes
advantage of the fact that E[S(X)]=Adot(θ), which is much faster than computing first moment
and correlations as currently done.

# Examples
```julia-repl
julia> n=[1,1,1]; θ=[1.0,-1.0,1.0,0.0];
julia> ∇logZ = CpelAsm.get_grad_logZ(n,θ);
julia> CpelAsm.comp_mmlv(n,∇logZ)
3-element Array{Float64,1}:
 0.88079708
 0.11920292
 0.88079708
```
"""
function comp_mmlv(n::Vector{Int64},∇logZ::Vector{Float64})::Vector{Float64}

    # Return
    return abs.(round.(0.5.*(1.0./n.*∇logZ[1:length(n)].+1.0);digits=8))

end # end comp_mmlv
"""
    `comp_corr(EX,EXX)`

Function that returns the correlation vector between consecutive CpG sites given the vector of
first order moments E[X] and second order moments E[XX].

# Examples
```julia-repl
julia> CpelAsm.comp_corr(CpelAsm.comp_ex([4],[0.0],0.0),CpelAsm.comp_exx([4],[0.0],0.0))
3-element Array{Float64,1}:
 0.0
 0.0
 0.0
```
"""
function comp_corr(ex::Vector{Float64},exx::Vector{Float64})::Vector{Float64}

    # Initialize vector of 1's
    ov = ones(Float64,length(ex)-1)

    # Return ρ
    return round.((exx.-ex[1:(end-1)].*ex[2:end]) ./ sqrt.((ov.-ex[1:(end-1)].^2) .*
    (ov.-ex[2:end].^2));digits=8)

end # end comp_corr
"""
    `comp_cov([N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that returns the covariance matrix of a methylation vector given the `[N1,...,NK]`,
`[α1,...,αK]`, `β`, E[X], and E[XX].

# Examples
```julia-repl
julia> n=[4]; α=[1.0]; β=1.0;
julia> ex = CpelAsm.comp_ex(n,α,β);
julia> exx = CpelAsm.comp_exx(n,α,β);
julia> CpelAsm.comp_cov(n,α,β,ex,exx)
4×4 Array{Float64,2}:
 0.419974     0.0          4.44089e-16  3.33067e-16
 0.0          0.419974     1.11022e-16  4.44089e-16
 4.44089e-16  1.11022e-16  0.419974     0.0
 3.33067e-16  4.44089e-16  0.0          0.419974
```
"""
function comp_cov(n::Vector{Int64},a::Vector{Float64},b::Float64,ex::Vector{Float64},
                  exx::Vector{Float64})::Array{Float64,2}

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
matrix Σ. The components of the vector with the largest magnitude highlight the CpG sites in the
region with the largest amount of variability.

# Examples
```julia-repl
julia> n=[4]; α=[1.0]; β=1.0;
julia> ex = CpelAsm.comp_ex(n,α,β);
julia> exx = CpelAsm.comp_exx(n,α,β);
julia> Σ = CpelAsm.comp_cov(n,α,β,ex,exx);
julia> CpelAsm.comp_evec(Σ)
4-element Array{Float64,1}:
 0.6867257391267225
 0.16854601514378278
 0.16854601514379314
 0.6867257391267322
```
"""
function comp_evec(cov::Array{Float64,2})::Vector{Float64}

    # Return evec
    return abs.(eigvecs(cov)[:,size(cov)[1]])

end # end comp_evec
###################################################################################################
# VERIFICATION FUNCTIONS (NOT USED)
###################################################################################################
"""
    `comp_nme_xcal(Z,[N1,...,NK],[α1,...,αK],β,EX,EXX)`

Function that computes normalized methylation entropy (NME) over the homozygous CpG sites. This is
done by assuming an Ising model for the allele-specific methylation state vector of size
`[N1,...,NK]`, parameters `[α1,...,αK]` and `β`. The homozygous part of the allele-specific vector
is determined by binary vector Z (i.e., via Hadamard product Z*X, where * is the Hadamard product).

# Examples
```julia-repl
julia> n=[10]; z=trues(sum(n)); α=[0.0]; β=0.0;
julia> CpelAsm.comp_nme_xcal(z,n,α,β)
1.0
```
"""
function comp_nme_xcal(z::BitArray{1},n::Vector{Int64},a::Vector{Float64},b::Float64)::Float64

    # Loop over 𝒳h
    h = 0.0
    w = zeros(Int64,sum(n))
    hom_ind = findall(y->y==true,z)
    xcal = generate_xcal(length(hom_ind))
    @inbounds for x in xcal
        w[hom_ind] =+ x
        lkhd = comp_lkhd(w,n,a,b)
        h += lkhd * log(lkhd)
        w[hom_ind] =- x
    end

    # Return
    return round(-h/(sum(z)*LOG2);digits=8)

end # end comp_nme_xcal
"""
    `comp_uc_xcal(Z1,Z2,N1,N2,theta1,theta2)`

Function that exactly computes the UC over the homozygous CpG sites. This is done by assuming an
Ising model for the allele-specific methylation state vector of size `[N1,...,NK]`, parameters
`[α1,...,αK]` and `β`. The homozygous part of the allele-specific vector is determined by binary
vector Z (i.e., via Hadamard product Z*X, where * is the Hadamard product).

# Examples
```julia-repl
julia> n1=[10]; n2=[10];
julia> z1=trues(sum(n1)); z2=trues(sum(n2));
julia> θ1=[2.0,0.0]; θ2=[-2.0,0.0];
julia> CpelAsm.comp_uc_xcal(z1,z2,n1,n2,θ1,θ2)
0.43482166257188215
```
"""
function comp_uc_xcal(z1::BitArray{1},z2::BitArray{1},n1::Vector{Int64},n2::Vector{Int64},
                      theta1::Vector{Float64},theta2::Vector{Float64})::Float64

    # Loop over 𝒳h
    num = 0.0
    den = 0.0
    n = sum(z1)
    w1 = zeros(Int64,sum(n1))
    w2 = zeros(Int64,sum(n2))
    hom_ind_1 = findall(y->y==true,z1)
    hom_ind_2 = findall(y->y==true,z2)
    @inbounds for i=1:(2^n)
        z = IntBitVec(i)[1:n]
        x = -ones(Int64,n)
        x[z] .= 1
        w1[hom_ind_1] =+ x
        w2[hom_ind_2] =+ x
        lkhd1 = comp_lkhd(w1,n1,theta1[1:(end-1)],theta1[end])
        lkhd2 = comp_lkhd(w2,n2,theta2[1:(end-1)],theta2[end])
        lkhd0 = 0.5*(lkhd1+lkhd2)
        den += lkhd0 * log(lkhd0)
        num += lkhd1 * log(lkhd1) + lkhd2 * log(lkhd2)
        w1[hom_ind_1] =- x
        w2[hom_ind_2] =- x
    end

    # Return
    return 1-0.5*num/den

end # end comp_uc_xcal
###################################################################################################
# COMPETING MODELS
###################################################################################################
"""
    `comp_lkhd_bin(XOBS,θ)`

Function that returns likelihood vector of XOBS, potentially partial, given parameter vector θ for
the non-parametric independent model (binomial).

"""
function comp_lkhd_bin(xobs::Array{Vector{Int64},1},θ::Vector{Float64})::Vector{Float64}

    # Return evec
    return [prod(θ.^(x.==1) .* (ones(length(xobs[1]))-θ).^(x.==-1)) for x in xobs]

end # end comp_lkhd_bin
"""
    `comp_lkhd_mult(XOBS,θ)`

Function that returns likelihood vector of XOBS, assuming full data, given parameter vector θ for
the non-parametric dependent model (multinomial).

"""
function comp_lkhd_mult(xobs::Array{Vector{Int64},1},θ::Vector{Float64})::Vector{Float64}

    # Generate 𝒳
    xcal = generate_xcal(length(xobs[1]))

    # Return evec
    return θ[[findfirst(y->y==x,xcal) for x in xobs]]

end # end comp_lkhd_mult
