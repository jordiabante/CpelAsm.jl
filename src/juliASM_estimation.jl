"""
    create_Ux(alpha,beta)

Function that generates a function to compute potential energy for
a region with N CpG sites and with parameters α and β.

# Examples
```julia-repl
julia> Ux_fun = create_Ux(N,alpha,beta)
```
"""
function create_Ux(a::Float64, b::Float64)
    function Ux(x::Array{Int64,1})::Float64
        -a * sum(x) - b * sum(x[2:end] .* x[1:(end-1)])
    end
    return Ux
end # end create_Ux
"""
comp_full_stats(xobs)

Compute vector S of sufficient statistics under model with N CpG sites.
This function is only valid for fully observed data since we are summing
over to obtain overall statistic S=[S1,S2] where
S1 = ∑ T1(x); S2 = ∑ T2(x)
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
comp_Z(N,alpha,beta)

Compute partition function of a model with N CpG cites and with parameters
α and β. This is computed using the closed-form expression provided in
REFERENCE. The transfer-matrix method is used to obtain a computationally
efficient expression without the need of recursive expressions.

# Examples
```julia-repl
julia>JuliASM.comp_Z(4,1.0,1.0)
1149.715905067999
```
"""
function comp_Z(N::Int64, a::Float64, b::Float64)::Float64
    # Eigenvalues
    aux1 = exp(b) * cosh(a)
    aux2 = sqrt(1 + exp(4*b)*sinh(a)^2)
    LN = [(aux1-exp(-b)*aux2)^(N-1) 0.0; 0.0 (aux1+exp(-b)*aux2)^(N-1)]

    # Matrix S
    aux1 = -exp(2*b) * sinh(a)
    S = [aux1-aux2 aux1+aux2; 1.0 1.0]

    # Matrix Sinv
    aux3 = 1.0 / aux2
    Sinv = [-0.5*aux3 0.5+0.5*aux1*aux3;0.5*aux3 0.5-0.5*aux1*aux3]

    # Return a Z>0
    u = [exp(-a/2.0); exp(a/2.0)]
    return max(1e-100, u' * S * LN * Sinv * u)
end
"""
comp_scal_fac(k,alpha,beta,alphap)

Compute scaling factor in a model with k unobserved CpG cites with parameters
α and β. This is computed using the closed-form expression provided in
REFERENCE. The transfer-matrix method is used to obtain a computationally
efficient expression without the need of recursive expressions.

# Examples
```julia-repl
julia>JuliASM.comp_scal_fac(1,1.0,1.0,1.0)
4.7048192304864935
```
"""
function comp_scal_fac(k::Int64, a::Float64, b::Float64, ap::Float64)::Float64
    # Return one if k=0
    k==0 && return 1.0

    # Eigenvalues
    aux1 = exp(b) * cosh(a)
    aux2 = sqrt(1 + exp(4*b)*sinh(a)^2)
    LN = [(aux1-exp(-b)*aux2)^(k-1) 0.0; 0.0 (aux1+exp(-b)*aux2)^(k-1)]

    # Matrix S
    aux1 = -exp(2*b) * sinh(a)
    S = [aux1-aux2 aux1+aux2; 1.0 1.0]

    # Matrix Sinv
    aux3 = 1.0 / aux2
    Sinv = [-0.5*aux3 0.5+0.5*aux1*aux3;0.5*aux3 0.5-0.5*aux1*aux3]

    # Return scaling factor
    vs = [exp(-ap); exp(ap)]
    ve = [exp(-a/2.0); exp(a/2.0)]
    return max(1e-100, vs' * S * LN * Sinv * ve)
end
"""
comp_lkhd(x,alpha,beta)

Compute likelihood of a partial/full observation x that can be missing
values on both sides (given by 0's). Returns -1 if x is not a valid
observation.

# Examples
```julia-repl
julia> JuliASM.comp_lkhd([0;1;-1;1;0],2.0,2.0)
6.144020836828058e-6
```

"""
function comp_lkhd(x::Array{Int64,1}, a::Float64, b::Float64)::Float64
    # Avoid the rest if N=1
    length(x)==1 && return 0.5 * exp(x[1]*a) / cosh(a)

    # Find kl and kr
    ind_dat = findall(!isequal(0), x)
    kl = ind_dat[1]-1
    kr = length(x)-ind_dat[end]

    # Get partition function and energy function
    Ux = create_Ux(a, b)
    Z = comp_Z(length(x), a, b)

    # Get scaling factors due to missing data
    sf = comp_scal_fac(kl, a, b, x[ind_dat[1]]*b+a/2.0) *
         comp_scal_fac(kr, a, b, x[ind_dat[end]]*b+a/2.0)

    # Return energy function evaluated at x properly scaled
    return exp(-Ux(x[ind_dat])) * sf / Z
end
"""
create_Llkhd(x,alpha,beta)

Create function to compute the minus log-likelihood function for a region with N
CpG sites given the M partial observations.

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
            ind_dat = findall(!isequal(0), x)
            kl = ind_dat[1]-1
            kr = length(x)-ind_dat[end]

            # Check if observation are contiguous
            all(ind_dat[2:end]-ind_dat[1:(length(ind_dat)-1)].==1) || continue

            # Add scaling factors due to missing data
            kl>0 && (aux+=log(comp_scal_fac(kl, a, b, x[ind_dat[1]]*b+a/2.0)))
            kr>0 && (aux+=log(comp_scal_fac(kr, a, b, x[ind_dat[end]]*b+a/2.0)))

            # Add log of it to Llkhd
            aux += -Ux(x[ind_dat])
        end
        # Add as many logZ as samples we have
        -aux + length(xobs) * logZ
    end
    return Llkhd_fun
end # end create_Llkhd
"""
est_alpha(xobs)

Estimate parameter α for the N=1 case.

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
est_eta(xobs)

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
mle_mult(xobs)

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a multinomial model
and full observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_mult(xobs)
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
mle_bin_semi(xobs)

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a binomial model that
assumes a different probability of methylation at each CpG site, as well as
fully observed vectors. This results in a semi-parametric model that grows
in size proportionally to the number of CpG sites.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_bin_semi(xobs)
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
mle_bin_param(xobs)

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a binomial model
that assumes that each CpG site has the same probability of being methylated,
as well as fully observed vectors. This results in a parametric model that does
not grow with the number of CpG sites considered.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_bin_param(xobs)
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
bifurcate(xpool,sel)

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
perm_test(xobs1,xobs2,tobs)

# Examples
```julia-repl
julia> ??
```
"""
function perm_test(xobs1::Array{Array{Int64,1},1},xobs2::Array{Array{Int64,1},1},
                   tobs::Float64,n::Int64)::Float64
    # Initialize
    better = worse = 0.0
    xpool = vcat(xobs1, xobs2)

    # Loop over combinations and compute statistic
    i = 1
    for subset in combinations(1:length(xpool), length(xobs2))
      # Bifurcate
      test, control = bifurcate(xpool, subset)

      # Estimate parameters
      n==1 ? eta1=est_alpha(control) : eta1=est_eta(control)
      n==1 ? eta2=est_alpha(test) : eta2=est_eta(test)

      # Compute MI for partition and compared to that observed
      comp_mi(n,eta1,eta2)>tobs ? better += 1.0 : worse += 1.0

      # If enough permutations leave
      i<100 ? i+=1 : break
    end

    # Return p-value
    return better/(worse+better)
end # end perm_test
