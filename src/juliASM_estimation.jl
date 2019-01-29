"""
    `create_Ux(α,β)`

Function that creates a function to compute potential energy for a region with
N CpG sites and with parameters α and β.

# Examples
```julia-repl
julia> Ux_fun = create_Ux(1.0,-1.0)
```
"""
function create_Ux(a::Float64, b::Float64)
    function Ux(x::Array{Int64,1})::Float64
        -a * sum(x) - b * sum(x[2:end] .* x[1:(end-1)])
    end
    return Ux
end # end create_Ux
"""
    `comp_full_stats(XOBS)`

Compute vector S of sufficient statistics under model with N CpG sites.
This function is only valid for fully observed data since we are summing
over to obtain overall statistic `S=[S1,S2]`` where
`S1 = ∑ T1(x)`; `S2 = ∑ T2(x)`
from all complete observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_full_data(10,4);
julia> Svec=JuliASM.comp_full_stats(xobs)
2-element Array{Int64,1}:
-4
-6
```
"""
function comp_full_stats(xobs::Array{Array{Int64,1},1})::Array{Int64,1}
    M = size(xobs)[1]
    s1 = sum(sum(xobs))
    s2 = sum([sum(xobs[i][2:end] .* xobs[i][1:(end-1)]) for i in 1:M])
    return [s1,s2]
end # comp_full_stats
"""
    `comp_Z(N,α,β)`

Compute partition function of a model with N CpG cites and with parameters
α and β. This expression is further simplified using the sum of rank-1 matrices
to avoid matrix multiplications.

# Examples
```julia-repl
julia>JuliASM.comp_Z(4,1.0,1.0)
1149.715905067998
```
"""
function comp_Z(N::Int64, a::Float64, b::Float64)::Float64
    # Eigenvalues
    aux1 = exp(b) * cosh(a)
    aux2 = sqrt(1.0 + exp(4.0*b)*sinh(a)^2)
    lambda1N = (aux1-exp(-b)*aux2)^(N-1)
    lambda2N = (aux1+exp(-b)*aux2)^(N-1)

    # Eigenvectors
    aux1 = -exp(2*b) * sinh(a)
    e1 = [aux1-aux2; 1.0]
    e1 /= sqrt(e1'*e1)
    e2 = [aux1+aux2; 1.0]
    e2 /= sqrt(e2'*e2)

    # Boundary conditions
    u = [exp(-a/2.0); exp(a/2.0)]

    # Return Z
    return max(1e-100,(u'*e1)^2*lambda1N + (u'*e2)^2*lambda2N)
end
"""
    `comp_scal_fac(K,α,β,αp1,αp2)`

Compute scaling factor in a model with K unobserved CpG cites with parameters
α and β. The transfer-matrix method is used to obtain a computationally
efficient expression without the need of recursive expressions. αp1 and αp2
are determined by the respective kind of boundary.

# Examples
```julia-repl
julia>JuliASM.comp_scal_fac(1,1.0,1.0,1.0,1.0)
7.524391382167261
```
"""
function comp_scal_fac(k::Int64, a::Float64, b::Float64, ap1::Float64, ap2::Float64)::Float64
    # Return one if k=0
    k==0 && return 1.0

    # Eigenvalues
    aux1 = exp(b) * cosh(a)
    aux2 = sqrt(1.0 + exp(4.0*b)*sinh(a)^2)
    lambda1k = (aux1-exp(-b)*aux2)^(k-1)
    lambda2k = (aux1+exp(-b)*aux2)^(k-1)

    # Eigenvectors
    aux1 = -exp(2*b) * sinh(a)
    e1 = [aux1-aux2; 1.0]
    e1 /= sqrt(e1'*e1)
    e2 = [aux1+aux2; 1.0]
    e2 /= sqrt(e2'*e2)

    # Boundary conditions
    vs = [exp(-ap1); exp(ap1)]
    ve = [exp(-ap2); exp(ap2)]

    # Return scaling factor
    return max(1e-100, vs'*e1*ve'*e1*lambda1k + vs'*e2*ve'*e2*lambda2k)
end
"""
    `comp_lkhd(X,α,β)`

Compute likelihood of a partial/full observation X that can be missing
values anywhere (given by 0's).

# Examples
```julia-repl
julia> JuliASM.comp_lkhd([0;1;-1;1;0],2.0,2.0)
6.144020836828058e-6
```

"""
function comp_lkhd(x::Array{Int64,1}, a::Float64, b::Float64)::Float64
    # Avoid the rest if N=1
    length(x)==1 && return 0.5 * exp(x[1]*a) / cosh(a)

    # Get partition function and energy function
    Ux = create_Ux(a, b)
    Z = comp_Z(length(x), a, b)

    # Find changes to/from 0 in x vector
    ind_zero_st = findall(isequal(-1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1
    ind_zero_end = findall(isequal(1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1

    # Determine whether it starts/finishes with 0
    x[1]==0 && pushfirst!(ind_zero_st,1)
    x[end]==0 && push!(ind_zero_end,length(x)+1)

    # Get scaling factors due to all missing data
    sf = 1.0
    for i in 1:length(ind_zero_st)
        k = ind_zero_end[i]-ind_zero_st[i]
        ind_zero_st[i]==1 ? ap1 = a/2.0 : ap1=x[ind_zero_st[i]-1]*b+a/2.0
        ind_zero_end[i]==length(x)+1 ? ap2 = a/2.0 : ap2=x[ind_zero_end[i]]*b+a/2.0
        sf *= comp_scal_fac(k, a, b, ap1, ap2)
    end

    # Return energy function properly scaled.
    # Note: when we have missing data the zeros don't contribute at all.
    return exp(-Ux(x)) * sf / Z
end
"""
    `create_Llkhd(XOBS)`

Create function to compute the minus log-likelihood function for a region with N
CpG sites given the M partial observations XOBS.

# Examples
```julia-repl
julia>xobs=JuliASM.gen_ising_full_data(10,4);
julia>LogLike=JuliASM.create_Llkhd(xobs)
```
"""
function create_Llkhd(xobs::Array{Array{Int64,1},1})
    function Llkhd_fun(eta::Array{Float64,1})::Float64
        # Get parameters
        aux = 0.0
        a = eta[1]
        b = eta[2]

        # Get energy function and partition function
        Ux = create_Ux(a, b)
        logZ = log(comp_Z(length(xobs[1]), a, b))

        # Contribution of each observation
        for x in xobs

            # Find changes to/from 0 in x vector
            ind_zero_st = findall(isequal(-1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1
            ind_zero_end = findall(isequal(1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1

            # Determine whether it starts/finishes with 0
            x[1]==0 && pushfirst!(ind_zero_st,1)
            x[end]==0 && push!(ind_zero_end,length(x)+1)

            # Get scaling factors due to all missing data
            for i in 1:length(ind_zero_st)
                k = ind_zero_end[i]-ind_zero_st[i]
                ind_zero_st[i]==1 ? ap1 = a/2.0 : ap1=x[ind_zero_st[i]-1]*b+a/2.0
                ind_zero_end[i]==length(x)+1 ? ap2 = a/2.0 : ap2=x[ind_zero_end[i]]*b+a/2.0
                aux += log(comp_scal_fac(k, a, b, ap1, ap2))
            end

            # Add log of it to Llkhd
            aux += -Ux(x)
        end

        # Return MINUS log-likelihood. Add as many logZ as samples we have
        -aux + length(xobs) * logZ
    end
    return Llkhd_fun
end # end create_Llkhd
"""
    `est_alpha(XOBS)`

Estimate parameter α for N=1 case.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_ising_full_data(100,1);
julia> JuliASM.est_alpha(xobs)
-0.020002667306849575
```
"""
function est_alpha(xobs::Array{Array{Int64,1},1})::Array{Float64,1}
    # Derivation of estimate:
        # log(p)+log(2) = α-log(cosh(α))
        # log(1-p)+log(2) = -α-log(cosh(α))
        # α = (log(p)-log(1-p))/2

    # Proportion of X=1
    phat = length(findall(x->x==[1],xobs)) / length(xobs)

    # Return estimate
    return [0.5*(log(phat)-log(1.0-phat)),0.0]
end # end est_alpha
"""
    `est_eta(XOBS)`

Estimate parameter vector η=[α, β] based on full or partial observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_ising_full_data(100,4);
julia> JuliASM.est_eta(xobs)
2-element Array{Float64,1}:
-0.021471185237893084
-0.04713465466621405
```
"""
function est_eta(xobs::Array{Array{Int64,1},1})::Array{Float64,1}
    # Handle missing data in optimization
    etahat = [NaN, NaN]
    if !(length(findall(x->0 in x, xobs))>0)
        # Should have a global minimum. Note: NelderMead can't be boxed
        eta0 = [0.0, 0.0]
        L = create_Llkhd(xobs)
        optim = optimize(L,eta0,NelderMead(),Optim.Options(iterations=100);
                         autodiff=:forward)
        etahat = optim.minimizer
    else
        # Multiple minima: try different initializations w/ NelderMead
        aux = Inf
        L = create_Llkhd(xobs)
        etas = [[.0,.0],[-.5,-.5],[.5,-.5],[-.5,.5],[.5,.5]]
        for eta0 in etas
            optim = optimize(L,eta0,NelderMead(),Optim.Options(iterations=100);
                             autodiff=:forward)
            if optim.minimum < aux
                aux = optim.minimum
                etahat = optim.minimizer
            end
        end
        # Boxed Simulated Annealing (SAMI)
        # lower = [-2.0, -2.0]
        # upper = [2.0, 2.0]
        # optim = optimize(L, lower, upper, eta0, SAMIN(rt=0.5; f_tol=1e-10,
        #     x_tol=1e-2), Optim.Options(iterations=10^6))
        # etahat = optim.minimizer
        # # BlackBoxOptim: global optimization package.
        # L = create_Llkhd(xobs)
        # optim = bboptimize(L; SearchRange = (-2.0, 2.0), NumDimensions = 2,
        #     MaxTime = 5.0, TraceMode=:silent)
        # # Return estimate vector
        # return best_candidate(optim)
    end
    # Return estimate
    return etahat
end # end est_eta
"""
    `mle_mult(XOBS)`

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a multinomial model
and full observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_mult(xobs)
16-element Array{Float64,1}:
 0.07
 0.08
 0.06
 0.06
 0.08
 0.08
 0.04
 0.05
 0.04
 0.07
 0.07
 0.02
 0.04
 0.12
 0.07
 0.05
```
"""
function mle_mult(xobs::Array{Array{Int64,1},1})::Array{Float64,1}

    # Get xcal
    N = size(xobs[1])[1]
    xcal = generate_xcal(N)

    # Get empirical estimates
    phat = zeros(Float64,2^N)
    for x in xobs
        phat[findfirst(y-> x==y, xcal)] += 1
    end
    phat = phat ./ sum(phat)

    # Return phat
    return phat
end # end mle_mult
"""
    `mle_bin_semi(XOBS)`

Estimate parameter vector `p=[p_1,...,p_{2^N}]` based on a binomial model that
assumes a different probability of methylation at each CpG site, as well as
fully observed vectors. This results in a semi-parametric model that grows
in size proportionally to the number of CpG sites.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_bin_semi(xobs)
16-element Array{Float64,1}:
 0.06662343999999999
 0.07512856000000001
 0.04824455999999999
 0.05440344
 0.07512856000000001
 0.08471944000000002
 0.054403440000000004
 0.06134856
 0.061498559999999994
 0.06934944000000001
 0.044533439999999994
 0.050218559999999995
 0.06934944000000001
 0.07820256000000002
 0.05021856
 0.056629439999999996
```
"""
function mle_bin_semi(xobs::Array{Array{Int64,1},1})::Array{Float64,1}

    # Get xcal
    N = size(xobs[1])[1]
    M = length(xobs)
    xcal = generate_xcal(N)

    # Get empirical estimate
    θ = zeros(Float64,N)
    [θ[i] = length(findall(x->x[i]==1,xobs)) / M for i in 1:N]

    # Compute probabilities of each pattern
    i = 1
    phat = zeros(Float64,2^N)
    for x in xcal
        phat[i] =prod([θ[j]^isequal(x[j],1) * (1-θ[j])^(1-isequal(x[j],1))
            for j in 1:N])
        i += 1
    end

    # Return phat
    return phat #/ sum(phat)
end # end mle_bin_semi
"""
    `mle_bin_param(XOBS)`

Estimate parameter vector `p=[p_1,...,p_{2^N}]` based on a binomial model
that assumes that each CpG site has the same probability of being methylated,
as well as fully observed vectors. This results in a parametric model that does
not grow with the number of CpG sites considered.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_bin_param(xobs)
16-element Array{Float64,1}:
 0.06765201
 0.06499899
 0.06499899
 0.06245000999999999
 0.06499899
 0.06245000999999999
 0.06245000999999999
 0.06000099
 0.06499899
 0.06245000999999999
 0.06245000999999999
 0.06000099
 0.06245000999999999
 0.06000099
 0.06000099
 0.05764800999999999
```

"""
function mle_bin_param(xobs::Array{Array{Int64,1},1})::Array{Float64,1}

    # Get xcal
    N = size(xobs[1])[1]
    M = length(xobs)
    xcal = generate_xcal(N)

    # Get empirical estimates θ=[θ_1,...,θ_N]
    θ = sum([length(findall(x->x==1,x)) / N for x in xobs]) / M

    # Compute probabilities of each pattern
    i = 1
    phat = zeros(Float64,2^N)
    for x in xcal
        k = length(findall(y->y==1,x))
        phat[i] = θ^(k) * (1-θ)^(N-k)
        i+=1
    end

    # Return phat
    return phat #/ sum(phat)
end # end mle_bin_param
"""
    `bifurcate(xpool,sel)`

Function that divides xpool into two mutually exclusive sets based on sel.

# Examples
```julia-repl
julia> JuliASM.bifurcate([[1,1],[1,-1],[-1,1]], [1,2])
(Array{Int64,1}[[1, 1], [1, -1]], Array{Int64,1}[[-1, 1]])
```
"""
function bifurcate(xpool::Array{Array{Int64,1},1}, sel::Vector{T}) where T <: Integer
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
